from sac_algo.utils import get_parsed_metrics
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams['mathtext.fontset'] = 'stix'
mpl.rc('font', family='times new roman')
mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['font.size'] = 10

data_ = get_parsed_metrics("train_log/SAC_test_8")

_, ax = plt.subplots(1, 3, figsize=(7.1, 1.95), constrained_layout=True)

ax[0].plot(
    data_["rollout/ep_rew_mean"][:,0],
    data_["rollout/ep_rew_mean"][:,1]
)

ax[0].set_xlabel("timesteps\n(a)")
ax[0].set_ylabel("mean episodic reward")

ax[0].grid(alpha=0.3, linestyle="-.")

ax[0].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)


ax[1].plot(
    data_["train/entropy"][:,0],
    data_["train/entropy"][:,1]
)

ax[1].set_xlabel("timesteps\n(b)")
ax[1].set_ylabel("entropy")

ax[1].grid(alpha=0.3, linestyle="-.")

ax[1].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)


ax[2].plot(
    data_["train/std1"][:,0],
    data_["train/std1"][:,1],
    label="$i_{invd}$"
)

ax[2].plot(
    data_["train/std2"][:,0],
    data_["train/std2"][:,1],
    label="$i_{invq}$"
)

ax[2].set_xlabel("timesteps\n(c)")
ax[2].set_ylabel("$\sigma$")

ax[2].grid(alpha=0.3, linestyle="-.")
ax[2].legend(loc=0, labelspacing=0.15, borderpad=0.15)
ax[2].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)

plt.savefig("Results/results_metrics.svg", dpi=600)
plt.show()