import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

newparams = {'axes.grid': False, 'lines.linewidth': 1,
             'font.size': 13, 'mathtext.fontset': 'stix'}

plt.rcParams.update(newparams)  # Updates the parameters set above.

def plot_taus_temp():
    taus_50 = np.genfromtxt("tau_extended_50.csv", delimiter=',')[1:]
    temps_50 = np.genfromtxt("temps_extended_50.csv", delimiter=',')[1:]
    taus_100 = np.genfromtxt("tau_extended_100.csv", delimiter=',')
    temps_100 = np.genfromtxt("temps_extended_100.csv", delimiter=',')
    taus_10 = np.genfromtxt("tau_extended_10.csv", delimiter=',')
    temps_10 = np.genfromtxt("temps_extended_10.csv", delimiter=',')
    taus_20 = np.genfromtxt("tau_extended_20.csv", delimiter=',')
    temps_20 = np.genfromtxt("temps_extended_20.csv", delimiter=',')
    T_2 = np.linspace(0, 2.23, 10000)
    m_analytical = (1 - (np.sinh(2 / T_2)) ** (-4)) ** (1 / 8)*2
    T_c = 2.2691
    t = np.linspace(0, 1)
    plt.plot(temps_10, taus_10, label=r'$N=10$')
    plt.plot(temps_10, taus_10, 'o')
    plt.plot(temps_20, taus_20, label=r'$N=20$')
    plt.plot(temps_20, taus_20, 'o')
    plt.plot(temps_50, taus_50,label=r'$N=50$')
    plt.plot(temps_50, taus_50, 'o')
    plt.plot(temps_100, taus_100, label=r'$N=100$')
    plt.plot(temps_100, taus_100, 'o')

    #plt.plot(T_2, m_analytical, label="Analytical")
    plt.axvline(T_c, alpha=0.2)
    plt.xlabel(r'$T$')
    #plt.title("N=100")
    plt.legend()
    plt.ylabel(r'$\tau / J$')
    #plt.show()
    #plt.savefig("figures/tau_T_multiple_N_extended.pdf", dpi=200)

def plot_convergence_mag():

    mag_150 = np.genfromtxt("conv_150.csv", delimiter=',')[:6000]
    mag = np.genfromtxt("conv_100.csv", delimiter=',')[:len(mag_150)]
    mag_70 = np.genfromtxt("conv_70.csv", delimiter=',')[:len(mag_150)]
    #mag_200 = np.genfromtxt("conv_200.csv", delimiter=',')
    sweeps = np.arange(0, len(mag))
    sweeps_ = np.arange(0, len(mag_150))
    plt.plot(sweeps, 2*mag_70/(70*70), label="$N=70$")
    plt.plot(sweeps, 2 * mag / (100 * 100), label="$N=100$")
    plt.plot(sweeps_, 2*mag_150/(150*150), label="$N=150$")
    #plt.plot(sweeps, 2*mag_200/(200*200), label="$N=200$")
    plt.xlabel('t')
    plt.ylabel(r'$\sigma$')
    plt.title("$T=1.9 K$")
    plt.legend()
    #plt.savefig("figures/conv_test_N_change.pdf", dpi=200)
    plt.show()


def plot_Ntau_N():

    #tauN_extended = np.genfromtxt("tauN_extended.csv", delimiter=',')
    #N_extended  = np.arange(0, len(tauN_extended))
    #N_extended = np.genfromtxt("N_extended.csv", delimiter=',')
    tauN = np.genfromtxt("tauN_ext.csv", delimiter=',')
    N= np.genfromtxt("N_ext.csv", delimiter=',')
    #plt.plot(N, intercept + slope*N, label="Least Squares")
    #print(slope)
    log_tau = np.log(tauN/N)
    log_N = np.log(N)
    #plt.plot(log_N, log_tau)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_N, log_tau)
    print(slope)

    plt.loglog(N, tauN/N, 'o', label="Ordinary")
    #plt.plot(N, tauN, 'o')
    #plt.plot(N_extended, tauN_extended, 'o', label="Extended")
    plt.xlabel(r"$N$")
    plt.ylabel(r"$\tau$")
    plt.title(r"$T = T_c$")
    plt.legend()
    #plt.savefig("figures/N_tauN_ord.pdf", dpi=200)
    plt.show()

def plot_lattice():
    """
    Visualizes the lattice, either before or after MC sweeps
    """
    fig = plt.figure(figsize=(6, 6))
    before = np.genfromtxt("lattice_before.csv", delimiter=',')
    after = np.genfromtxt("lattice_after.csv", delimiter=',')
    plt.imshow(before, cmap='inferno', interpolation='none', )

    plt.xticks([])
    plt.title("Initial state, $N=100, T = 1.9 K$")
    plt.yticks([])
    #plt.show()
    #plt.savefig("figures/lattice_before.pdf", dpi=500)


def plot_tau_o_N():
    tau_1 = np.genfromtxt("tau_1.csv", delimiter=',')
    N_1 = np.genfromtxt("N_1.csv", delimiter=',')
    t_1 = 0.56
    tau_15 = np.genfromtxt("tau_15.csv", delimiter=',')
    N_15 = np.genfromtxt("N_15.csv", delimiter=',')
    t_15 = 0.34
    tau_18 = np.genfromtxt("tau_18.csv", delimiter=',')
    N_18 = np.genfromtxt("N_18.csv", delimiter=',')
    t_18 = 0.2
    tau_21 = np.genfromtxt("tau_21.csv", delimiter=',')
    N_21 = np.genfromtxt("N_21.csv", delimiter=',')
    t_21 = 0.07
    tau_05 = np.genfromtxt("tau_05.csv", delimiter=',')
    N_05 = np.genfromtxt("N_05.csv", delimiter=',')
    t_05 = 0.78

    tau_0 = 3.99
    #tau_temp_2 = np.genfromtxt("tau_temp_2.csv", delimiter=',')
    #int_N_t_2 = np.genfromtxt("inv_N_t_2.csv", delimiter=',')
    plt.plot(1 / (N_05 * t_05), tau_05 / t_05, 'o', label=r"$t=0.78$")
    plt.plot(1/(N_1*t_1), tau_1/t_1, 'o', label=r"$t=0.55$")
    plt.plot(1/(N_15*t_15), tau_15/t_15, 'o', label=r"$t=0.40$")
    plt.plot(1 / (N_18 * t_18), tau_18 / t_18, 'o', label=r"$t=0.20$")
    #plt.plot(1 / (N_21 * t_21), tau_21/ t_21, 'o', label=r"$t=0.07$")
    plt.ylabel(r"$\tau/t$")
    plt.xlabel(r"$1/Nt$")
    plt.legend(loc=4)
    plt.show()
    #plt.savefig("figures/tau_not.pdf", dpi=500)

def plot_other():
    """
    Plots task e), for multiple values of t. Note; the files has to be generated first.
    """
    tau_1 = np.genfromtxt("tau_1.csv", delimiter=',')
    N_1 = np.genfromtxt("N_1.csv", delimiter=',')
    t_1 = 0.56
    tau_15 = np.genfromtxt("tau_15.csv", delimiter=',')
    N_15 = np.genfromtxt("N_15.csv", delimiter=',')
    t_15 = 0.34
    tau_18 = np.genfromtxt("tau_18.csv", delimiter=',')
    N_18 = np.genfromtxt("N_18.csv", delimiter=',')
    t_18 = 0.2
    tau_05 = np.genfromtxt("tau_05.csv", delimiter=',')
    N_05 = np.genfromtxt("N_05.csv", delimiter=',')
    t_05 = 0.78
    tau_223 = np.genfromtxt("tau_223.csv", delimiter=',')
    N_223 = np.genfromtxt("N_223.csv", delimiter=',')
    t_223 = 0.02
    tau_01 = np.genfromtxt("tau_01.csv", delimiter=',')
    N_01 = np.genfromtxt("N_01.csv", delimiter=',')
    t_01 = 0.96

    # Analytical value of tau, for comparison
    tau_0 = 3.99

    plt.plot((N_01 * t_01), tau_01 / t_01, 'o', label=r"$t=0.96$")
    plt.plot((N_05 * t_05), tau_05 / t_05, 'o', label=r"$t=0.78$")
    plt.plot((N_1 * t_1), tau_1 / t_1, 'o', label=r"$t=0.55$")
    plt.plot((N_15 * t_15), tau_15 / t_15, 'o', label=r"$t=0.40$")
    plt.plot((N_18 * t_18), tau_18 / t_18, 'o', label=r"$t=0.20$")
    plt.plot((N_223 * t_223), tau_223 / t_223, 'o', label=r"$t=0.02$")
    plt.axhline(tau_0)
    plt.ylabel(r"$\tau/t$")
    plt.xlabel(r"$Nt$")
    plt.legend()
    plt.show()
    # Computes the mean, to get  tau_0.
    print(np.mean(tau_15 / t_15))
    #plt.savefig("figures/tau_not_found.pdf", dpi=500)

plot_other()
#plot_tau_o_N()
#plot_lattice()
#plot_Ntau_N()
#plot_taus_temp()
#plot_convergence_mag()
