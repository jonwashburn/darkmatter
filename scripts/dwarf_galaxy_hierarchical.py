#!/usr/bin/env python3
"""
Hierarchical Bayesian analysis of dwarf galaxy kinematics using the
bandwidth-limited gravity framework. This approach provides robust parameter
estimation by defining priors and using MCMC to explore the full posterior.
"""
import numpy as np
import pandas as pd
import emcee
import pickle
from pathlib import Path
import matplotlib.pyplot as plt

# --- Model & Data ---

def load_dwarf_data():
    """Load and filter dwarf galaxy data."""
    data_path = Path("../data/dwarfs/dwarf_all.csv")
    df = pd.read_csv(data_path)
    required = ['vlos_sigma', 'rhalf', 'distance', 'mass_stellar']
    df = df.dropna(subset=required)
    df = df[(df['vlos_sigma'] > 0) & (df['rhalf'] > 0) & 
            (df['distance'] > 0) & (df['mass_stellar'] > 0)]
    df['vlos_sigma_err'] = df['vlos_sigma_em'].fillna(df['vlos_sigma'] * 0.1)
    return df

def predict_sigma(theta, galaxy_data):
    """Predicts sigma_vlos from model parameters."""
    lambda_bw, alpha_urgency, xi_base, log_sigma_intr = theta
    
    # Unpack data for clarity
    m_stellar = galaxy_data['mass_stellar']
    r_half_arcsec = galaxy_data['rhalf']
    distance_kpc = galaxy_data['distance']

    # Physical constants and conversions
    G_kpc = 4.302e-6
    c_light = 3e5
    tau_0 = 1e8
    
    # Calculate physical radius
    r_half_kpc = r_half_arcsec * distance_kpc / 206265.0
    
    # Total mass including the Branching Factor (fixed from first principles)
    BRANCHING_FACTOR = 828.4
    m_total = m_stellar * BRANCHING_FACTOR
    
    # Newtonian expectation
    sigma_newton = np.sqrt(G_kpc * m_total / (5.0 * r_half_kpc))
    
    # Bandwidth correction
    omega_dyn = sigma_newton / np.maximum(r_half_kpc, 1e-3)
    n_complexity = xi_base * (1.0 + 5.0 / np.maximum(r_half_kpc, 0.01))
    B_max = lambda_bw * c_light / (4 * np.pi * G_kpc)
    tau_update = n_complexity / (alpha_urgency * B_max)
    delay_factor = 1.0 / (1.0 + tau_update * omega_dyn / tau_0)
    
    return sigma_newton * delay_factor

# --- MCMC Framework ---

def log_prior(theta):
    """Log-prior for the model parameters."""
    lambda_bw, alpha_urgency, xi_base, log_sigma_intr = theta
    
    # Priors based on first principles & SPARC fits
    if not (0.01 < lambda_bw < 0.5): return -np.inf
    if not (0.1 < alpha_urgency < 5.0): return -np.inf
    if not (0.1 < xi_base < 10.0): return -np.inf
    if not (-3 < log_sigma_intr < 2): return -np.inf # Intrinsic scatter
        
    return 0.0

def log_likelihood(theta, data):
    """Log-likelihood of the data given the model."""
    sigma_pred = predict_sigma(theta, data)
    
    # Unpack observed data
    sigma_obs = data['vlos_sigma']
    sigma_err = data['vlos_sigma_err']
    
    # Total variance is measurement error + intrinsic scatter
    log_sigma_intr = theta[-1]
    sigma2 = sigma_err**2 + (10**log_sigma_intr)**2
    
    return -0.5 * np.sum((sigma_obs - sigma_pred)**2 / sigma2 + np.log(sigma2))

def log_probability(theta, data):
    """Full log-probability (prior + likelihood)."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, data)

# --- Main Execution ---

def main():
    """Main execution block."""
    # Load data
    data = load_dwarf_data()
    print(f"Loaded {len(data)} dwarf galaxies with required data.")

    # MCMC setup
    nwalkers = 32
    ndim = 4  # lambda_bw, alpha_urgency, xi_base, log_sigma_intr
    nsteps = 5000
    burnin = 1500

    # Initial guess for the sampler
    initial_guess = np.array([0.12, 1.0, 1.5, 0.0]) # λ, α, ξ, log_σ_intr
    pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)

    # Set up the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data,))
    
    print("Running MCMC...")
    sampler.run_mcmc(pos, nsteps, progress=True)

    # Get results
    flat_samples = sampler.get_chain(discard=burnin, thin=15, flat=True)
    
    # Print summary
    print("\nMCMC Results (median and 1-sigma uncertainties):")
    labels = ["lambda_bw", "alpha_urgency", "xi_base", "log_sigma_intr"]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(f"{labels[i]:>15}: {mcmc[1]:.3f} [-{q[0]:.3f}, +{q[1]:.3f}]")

    # Save results
    with open('dwarf_hierarchical_results.pkl', 'wb') as f:
        pickle.dump({'samples': flat_samples, 'labels': labels}, f)
    print("\nSaved results to dwarf_hierarchical_results.pkl")

    # Plotting (optional, corner plot)
    try:
        import corner
        fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True)
        fig.savefig("dwarf_corner_plot.png")
        print("Saved corner plot to dwarf_corner_plot.png")
    except ImportError:
        print("\nInstall 'corner' for plotting: pip install corner")

if __name__ == "__main__":
    main() 