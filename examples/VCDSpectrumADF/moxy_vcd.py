from scm.plams import AMSJob, Molecule
from scm.input_classes import drivers, engines

# import matplotlib.pyplot as plt

# This is a very simple script for the calculation of the VCD
# spectrum of (R)-methyloxirane using ADF at the PBE/DZP level
# of theory. Note that we use a pre-optimized geometry that was
# obtained at the same level of theory as the VCD calculation.
# If not already available you should include the optimization step.

# To plot the spectrum you can uncomment all lines related to matplotlib


def main():
    """
    The main script
    """

    # Define molecule, here we use a file with a pre-optimized geometry
    # Note: the geometry has been optimized at the same level of theory
    #       that we are now using for the VCD spectrum!!!
    methyloxirane = Molecule("R-methyloxirane.xyz")

    # Define engine
    driver = drivers.AMS()
    driver.Task = "SinglePoint"
    driver.Properties.VCD = "Yes"
    driver.Engine = engines.ADF()

    # Define settings
    driver.Engine.Basis.Core = "None"
    driver.Engine.Basis.Type = "DZP"
    driver.Engine.XC.GGA = "PBE"

    # Run job
    job = AMSJob(molecule=methyloxirane, settings=driver, name="moxy_vcd")
    results = job.run()

    # Collect results, first we get the rotational strengths
    R = results.get_vcd_rotational_strength()
    print(f"VCD rotational strengths in (10-44 esu^2 cm^2):\n {R}")
    # Then we print the rotational strength convoluted spectrum
    spectrum = results.get_vcd_spectrum()
    print(f"VCD rotational strengths spectrum:\n {spectrum}")
    # Optionally plot the spectrum
    # x_freq, y_intens_raw = spectrum
    # plt.plot(x_freq, y_intens_raw)
    # plt.xlabel("Frequency (cm^-1)")
    # plt.title("VCD rotational strength spectrum (in 10^(-44) esu^2 cm^2)")
    # plt.xlim(500, 4000)
    # plt.show()


if __name__ == "__main__":
    main()
