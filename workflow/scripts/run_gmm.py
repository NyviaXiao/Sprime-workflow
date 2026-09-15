"""Fit and compare one-, two- and three-component Gaussian mixtures.

Each invocation analyzes one population and one configured Denisovan match-rate
column. Model selection and component parameters are emitted as tidy TSV files.
"""
import numpy as np
from scipy.stats import chi2, norm
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import read, write, rate

SELECTION = ['population', 'reference', 'n_segments', 'status', 'method', 'selected_components',
             'loglik_1', 'loglik_2', 'loglik_3', 'lr_1_vs_2', 'p_1_vs_2', 'adjusted_p_1_vs_2',
             'second_comparison', 'lr_second', 'p_second', 'adjusted_p_second', 'converged']
COMPONENTS = ['population', 'reference', 'component', 'mean', 'standard_deviation', 'weight']


def fit(x, k, cfg, seed):
    """Fit a k-component model and return it with total log-likelihood."""
    model = GaussianMixture(n_components=k, n_init=cfg['n_init'], random_state=seed,
                            max_iter=500, reg_covar=1e-6).fit(x.reshape(-1, 1))
    if not model.converged_:
        raise RuntimeError(f'GMM k={k} did not converge')
    return model, model.score(x.reshape(-1, 1))*len(x)


def compare(x, low, high, models, cfg, seed):
    """Compare nested component counts using the configured calibration."""
    # GaussianMixture.score returns mean log-likelihood per observation; fit()
    # converted it to total log-likelihood before this LRT statistic.
    lr = max(0., 2*(models[high][1]-models[low][1]))
    if cfg['method'] == 'legacy_chi_square':
        return lr, float(chi2.sf(lr, 3*(high-low)))
    # Parametric bootstrap under the fitted lower-component null. The +1
    # correction prevents a zero p-value with a finite number of replicates.
    # This calibrates the mixture comparison, not upstream biological filtering.
    rng = np.random.default_rng(seed)
    null = models[low][0]
    means, std = null.means_.ravel(), np.sqrt(null.covariances_.ravel())
    exceed = 0
    b = int(cfg['bootstrap_replicates'])
    if b < 1:
        raise ValueError('bootstrap_replicates must be positive')
    for i in range(b):
        component = rng.choice(low, size=len(x), p=null.weights_)
        simulated = rng.normal(means[component], std[component])
        _, ll0 = fit(simulated, low, cfg, seed+i+1)
        _, ll1 = fit(simulated, high, cfg, seed+i+1)
        exceed += max(0., 2*(ll1-ll0)) >= lr
    return lr, (exceed+1)/(b+1)


def analyze(s):
    """Select a GMM, write tables and render its density diagnostic plot."""
    cfg = s.config['gmm']
    pop, ref = s.wildcards.population, s.wildcards.ref
    values = [rate(r, ref) for r in read(s.input[0])]
    x = np.asarray([v for v in values if v is not None], dtype=float)
    if np.any(~np.isfinite(x)) or np.any((x < 0) | (x > 1)):
        raise ValueError('Match rates must be finite and within [0,1]')
    # Pre-fill optional statistics so insufficient-data rows have the identical
    # schema as successful fits and can be concatenated safely.
    result = dict.fromkeys(SELECTION, 'NA')
    result.update(population=pop, reference=ref, n_segments=len(x), method=cfg['method'])
    components = []
    fig, ax = plt.subplots(figsize=(6, 5))
    if len(x) < max(3, cfg['min_segments']) or len(np.unique(x)) < 3:
        result['status'] = 'insufficient_data'
        ax.text(.5, .5, f'Insufficient data (n={len(x)})', ha='center', transform=ax.transAxes)
    else:
        # Fit each model once and reuse it for likelihood tests and plotting.
        models = {k: fit(x, k, cfg, cfg['seed']) for k in (1, 2, 3)}
        if any(models[k+1][1] < models[k][1]-1e-6 for k in (1, 2)):
            raise RuntimeError('Higher-component fit has lower likelihood; increase n_init')
        # The correction family is declared by configuration rather than by
        # successful fits, so data-poor groups do not shrink the correction.
        family = len(s.config['populations'])*len(cfg['target_references'])
        lr12, p12 = compare(x, 1, 2, models, cfg, cfg['seed']+10000)
        # Preserve the original decision tree: a significant 1-vs-2 result
        # proceeds to 2-vs-3; otherwise the second comparison is 1-vs-3.
        low = 2 if p12*family < cfg['alpha'] else 1
        lr_next, p_next = compare(x, low, 3, models, cfg, cfg['seed']+20000)
        k = 3 if p_next*family*2 < cfg['alpha'] else low
        result.update(status='ok', selected_components=k, lr_1_vs_2=lr12, p_1_vs_2=p12,
                      adjusted_p_1_vs_2=min(1., p12*family), second_comparison=f'{low}_vs_3',
                      lr_second=lr_next, p_second=p_next, adjusted_p_second=min(1., p_next*family*2), converged=True)
        result.update({f'loglik_{k}': models[k][1] for k in models})
        model = models[k][0]
        grid, density = np.linspace(0, 1, 500), np.zeros(500)
        ax.hist(x, bins=25, density=True, alpha=.5, edgecolor='black', linewidth=.5)
        # Mixture labels are arbitrary. Sorting by mean makes component numbers
        # stable and interpretable between runs and populations.
        for i, idx in enumerate(np.argsort(model.means_.ravel()), 1):
            mean = model.means_.ravel()[idx]
            sd, weight = np.sqrt(model.covariances_.ravel()[idx]), model.weights_[idx]
            curve = weight*norm.pdf(grid, mean, sd)
            density += curve
            ax.plot(grid, curve, '--')
            components.append(dict(zip(COMPONENTS, [pop, ref, i, mean, sd, weight])))
        ax.plot(grid, density, color='black')
        ax.text(.02, .98, f'k={k}; n={len(x)}\n1 vs 2: p={p12:.3g}\n{low} vs 3: p={p_next:.3g}',
                va='top', transform=ax.transAxes)
    ax.set(xlim=(0, 1), xlabel=f'Match rate to {ref}', ylabel='Density', title=pop)
    fig.tight_layout()
    fig.savefig(s.output.plot, dpi=250)
    plt.close(fig)
    write(s.output.selection, SELECTION, [result])
    write(s.output.components, COMPONENTS, components)
