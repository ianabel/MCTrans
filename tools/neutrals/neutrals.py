import numpy as np

K = 8.617E-5 #[eV/K]
MP = 1.673E-27 #[kg]
ME = 9.11e-31 #[kg]
QE = 1.602e-19 #[C]
BohrRadius = 5.29177e-11 # m
FineStructureConstant = 7.2973525694e-3

species = {'proton': {'mass': 1, 'charge': 1},
    'hydrogen': {'mass': 1, 'charge': 0},
    'electron': {'mass': ME / MP, 'charge': 1}}

class Species:
    def __init__(self, name):
        self.name = name
        self.charge = species[self.name]['charge']
        self.mass = species[self.name]['mass']

class CrossSection:
    def __init__(self, sigma, energy, particle, target, reaction):
        self.sigma = sigma
        self.energy = energy
        self.particle = particle
        self.target = target
        self.reaction = reaction
        self.reducedMass = particle.mass * target.mass / (particle.mass + target.mass)

class RateCoefficient(CrossSection):
    def __init__(self, rateCoefficient, temperature, crossSection):
        super().__init__(crossSection.sigma, crossSection.energy, crossSection.particle, crossSection.target, crossSection.reaction)
        self.rateCoefficient = rateCoefficient
        self.temperature = temperature

    def reactionRate(self, particleDensity, targetDensity, plasmaVolume):
        self.reactionRate = particleDensity * targetDensity * self.rateCoefficient * plasmaVolume

class Neutrals:
    def __init__(self, MachNumber, temperature, energy):
        # Mach Number is defined in relation to the sound speed
        # Temperature is in eV
        self.MachNumber = MachNumber
        self.temperature = temperature
        self.energy = energy
        self.crossSections = np.array([])
        self.rateCoefficientsCold = np.array([])
        self.rateCoefficientsHot = np.array([])

        self.computeCrossSections(self.energy)
        self.computeRateCoefficients()


    def add_CrossSection(self, crossSection):
        self.crossSections = np.append(self.crossSections, crossSection)

    def add_RateCoefficient(self, rateCoefficient, type):
        if type == 'hot':
            self.rateCoefficientsHot = np.append(self.rateCoefficientsHot, rateCoefficient)
        elif type == 'cold':
            self.rateCoefficientsCold = np.append(self.rateCoefficientsCold, rateCoefficient)
        else:
            print('Please enter either hot or cold for the neutral type')
            return

    def computeCrossSections(self, energy):
        self.electronImpactIonizationCrossSection = self.electronImpactIonizationCrossSection(energy)
        self.protonImpactIonizationCrossSection = self.protonImpactIonizationCrossSection(energy)
        self.chargeExchangeCrossSection = self.chargeExchangeCrossSection(energy)
        self.chargeExchangeCrossSection1993 = self.chargeExchangeCrossSection1993(energy)
        self.radiativeRecombinationCrossSection = self.radiativeRecombinationCrossSection(energy)

        self.add_CrossSection([self.electronImpactIonizationCrossSection, self.protonImpactIonizationCrossSection, self.chargeExchangeCrossSection, self.chargeExchangeCrossSection1993, self.radiativeRecombinationCrossSection])

    def computeRateCoefficients(self):
        self.electronImpactIonizationRateCoefficientCold = self.rateCoefficientCold(self.electronImpactIonizationCrossSection)
        self.protonImpactIonizationRateCoefficientCold = self.rateCoefficientCold(self.protonImpactIonizationCrossSection)
        self.chargeExchangeRateCoefficientCold = self.rateCoefficientCold(self.chargeExchangeCrossSection)
        self.chargeExchangeRateCoefficientCold1993 = self.rateCoefficientCold(self.chargeExchangeCrossSection1993)
        self.radiativeRecombinationRateCoefficientCold = self.rateCoefficientCold(self.radiativeRecombinationCrossSection)

        self.electronImpactIonizationRateCoefficientHot = self.rateCoefficientHot(self.electronImpactIonizationCrossSection)
        self.protonImpactIonizationRateCoefficientHot = self.rateCoefficientHot(self.protonImpactIonizationCrossSection)
        self.chargeExchangeRateCoefficientHot = self.rateCoefficientHot(self.chargeExchangeCrossSection)
        self.chargeExchangeRateCoefficientHot1993 = self.rateCoefficientHot(self.chargeExchangeCrossSection1993)
        self.radiativeRecombinationRateCoefficientHot = self.rateCoefficientHot(self.radiativeRecombinationCrossSection)

        self.add_RateCoefficient([self.electronImpactIonizationRateCoefficientCold, self.protonImpactIonizationRateCoefficientCold, self.chargeExchangeRateCoefficientCold, self.chargeExchangeRateCoefficientCold1993, self.radiativeRecombinationRateCoefficientCold], 'cold')
        self.add_RateCoefficient([self.electronImpactIonizationRateCoefficientHot, self.protonImpactIonizationRateCoefficientHot, self.chargeExchangeRateCoefficientHot, self.chargeExchangeRateCoefficientHot1993, self.radiativeRecombinationRateCoefficientHot], 'hot')

    def electronImpactIonizationCrossSection(self, energy):
        # Contribution from ground state
        # Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
        # Equation 1.2.1
        # e + H(1s) --> e + H+ + e
        # Accuracy is 10% or better
        reaction = 'e + H(1s) --> e + H+ + e'
        particle = Species('electron')
        target = Species('hydrogen')

        # Minimum energy of cross section in eV
        ionizationEnergy = 13.6
        minimumEnergySigma = ionizationEnergy

        fittingParamA = 0.18450
        fittingParamB =  np.array([-0.032226, -0.034539, 1.4003, -2.8115, 2.2986])

        sigma = np.zeros(len(energy))
        for i, E in enumerate(energy):
            if (E < minimumEnergySigma):
                sigma[i] = 0
            else:
                sum = 0.0
                for n, param in enumerate(fittingParamB):
                    sum += param * np.power(1 - ionizationEnergy / E, n)

                sigma[i] = 1.0e-13 / (ionizationEnergy * E) * (fittingParamA * np.log(E / ionizationEnergy) + sum)

        return CrossSection(sigma, energy, particle, target, reaction)

    def radiativeRecombinationCrossSection(self, energy):
        # From https://iopscience-iop-org.proxy-um.researchport.umd.edu/article/10.1088/1402-4896/ab060a
        # Igor A Kotelnikov and Alexander I Milstein 2019 Phys. Scr. 94 055403
        # Equation 9
        # H+ + e --> H + hν
        reaction = 'H+ + e --> H + hν'
        particle = Species('electron')
        target = Species('proton')

        Z = 1
        IonizationEnergy = 13.59844 # eV
        J_Z = np.power( Z, 2 ) * IonizationEnergy
        eta = np.sqrt( J_Z / energy )
        sigma = np.power( 2, 8 ) * np.power( np.pi * BohrRadius, 2 ) / 3 * np.power( eta, 6 ) * np.exp( - 4 * eta * np.arctan( 1 / eta ) ) / ( ( 1 - np.exp( -2 * np.pi * eta ) ) * np.power( np.power( eta, 2 ) + 1, 2 ) ) * np.power( FineStructureConstant, 3 )

        sigma *= 1e4 # convert to cm^2
        return CrossSection(sigma, energy, particle, target, reaction)


    def protonImpactIonizationCrossSection(self, energy):
        # Contribution from ground state
    	# Janev 1993, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION, Volume 4
    	# Equation 2.2.1
    	# H+ + H(1s) --> H+ + H+ + e
    	# Accuracy is 30% or better
        reaction = 'H+ + H(1s) --> H+ + H+ + e'
        particle = Species('proton')
        target = Species('hydrogen')

        # Minimum energy of cross section in keV
        minimumEnergySigma = 0.2
        TiKEV = energy / 1000

        A1 = 12.899
        A2 = 61.897
        A3 = 9.2731e3
        A4 = 4.9749e-4
        A5 = 3.9890e-2
        A6 = -1.5900
        A7 = 3.1834
        A8 = -3.7154

        sigma = np.zeros(len(energy))
        for i, E in enumerate(TiKEV):
            if (E < minimumEnergySigma):
            	sigma[i] = 0
            else:
            	# Energy is in units of keV
            	sigma[i] = 1e-16 * A1 * (np.exp(-A2 / E) * np.log(1 + A3 * E) / E + A4 * np.exp(-A5 * E) / (np.power(E, A6) + A7 * np.power(E, A8)))

        return CrossSection(sigma, energy, particle, target, reaction)

    def chargeExchangeCrossSection(self, energy):
        # Contribution from ground -> ground state
        # Janev 1987 3.1.8
        # H+ + H(1s) --> H(1s) + H+
        reaction = 'H+ + H(1s) --> H(1s) + H+'
        particle = Species('proton')
        target = Species('hydrogen')

        # Minimum energy of cross section in eV
        minimumEnergySigma_1s = 0.1
        minimumEnergySigma_2p = 19.0
        minimumEnergySigma_2s = 0.1

        sigma = np.zeros(len(energy))
        for i, E in enumerate(energy):
            if (E < minimumEnergySigma_1s):
                sigma_1s = 0
            else:
                sigma_1s = 0.6937e-14 * np.power(1 - 0.155 * np.log10(E), 2) / (1 + 0.1112e-14 * np.power(E, 3.3))


            # Janev 1987 3.1.9
            # p + H(1s) --> H(2p) + p
            aSigma_2p = np.array([-2.197571949935e+01, -4.742502251260e+01, 3.628013140596e+01, -1.423003075866e+01, 3.273090240144e+00, -4.557928912260e-01, 3.773588347458e-02, -1.707904867106e-03, 3.251203344615e-05])
            # Janev 1987 3.1.10
            # p + H(1s) --> H(2s) + p
            aSigma_2s = np.array([-1.327325087764e+04, 1.317576614520e+04, -5.683932157858e+03, 1.386309780149e+03, -2.089794561307e+02, 1.992976245274e+01, -1.173800576157e+00, 3.902422810767e-02, -5.606240339932e-04])

            # Contribution from ground -> 2p orbital
            if (E < minimumEnergySigma_2p):
                sigma_2p = 0
            else:
                sigma_2p = self.evaluateJanevCrossSectionFit(aSigma_2p, E)


            # Contribution from ground -> 2s orbital
            if (E < minimumEnergySigma_2s):
                sigma_2s = 0
            else:
                sigma_2s = self.evaluateJanevCrossSectionFit(aSigma_2s, E)

            sigma[i] = sigma_1s + sigma_2s + sigma_2p

        return CrossSection(sigma, energy, particle, target, reaction)

    def chargeExchangeCrossSection1993(self, energy):
        # We need the energy in projectile form, while the input is in CoM
        # Factor of two because of hydrogen
        energy *= 2

        # Contribution from ground -> ground state
        # Janev 1993 2.3.1
        # H+ + H(1s) --> H + H+
        # Janev 1993 2.3.2
        # H+ + H(n) --> H + H+
        reaction = 'H+ + H --> H + H+ (1993)'
        particle = Species('proton')
        target = Species('hydrogen')

        # Minimum energy of cross section in eV
        minimumEnergySigma_n1 = 0.12
        minimumEnergySigma_n2 = 10
        minimumEnergySigma_n3 = 10

        sigma = np.zeros(len(energy))
        for i, E in enumerate(energy):
            if (E < minimumEnergySigma_n1):
                sigma_n1 = 0
            else:
                energyKEV = E / 1000
                sigma_n1 = 1e-16 * 3.2345 * np.log( 235.88 / energyKEV + 2.3713 ) / ( 1 + 0.038371 * energyKEV + 3.8068e-6 * np.power( energyKEV, 3.5 ) + 1.1832e-10 * np.power( energyKEV, 5.4 ) )

            # Contribution from n=2 -> ground
            if (E < minimumEnergySigma_n2):
                sigma_n2 = 0
            else:
                energyKEV = E / 1000
                n = 2
                Etilde = energyKEV * n**2
                sigma_n2 = 1e-16 * 0.92750 * np.log( 6.5040e3 / Etilde + 20.699 ) / ( 1 + 1.3405e-2 * Etilde + 3.0842e-6 * np.power( Etilde, 3.5 ) + 1.1832e-10 * np.power( Etilde, 5.4 ) )


            # Contribution from n=3 -> ground
            if (E < minimumEnergySigma_n3):
                sigma_n3 = 0
            else:
                energyKEV = E / 1000
                n = 3
                Etilde = energyKEV * n**2
                sigma_n3 = 1e-16 * 0.37271 * np.log( 2.7645e6 / Etilde + 1.4857e3 ) / ( 1 + 1.5720e-3 * Etilde + 3.0842e-6 * np.power( Etilde, 3.5 ) + 1.1832e-10 * np.power( Etilde, 5.4 ) )

            sigma[i] = sigma_n1 + sigma_n2 + sigma_n3

        return CrossSection(sigma, energy, particle, target, reaction)

    def evaluateJanevCrossSectionFit(self, PolynomialCoefficients, energy):
        N_JANEV_COEFFS = 9
        if (len(PolynomialCoefficients) != N_JANEV_COEFFS):
            print("Janev uses fixed order fits, there are not " + str(N_JANEV_COEFFS) + " numbers. Something is wrong.")
        sum = 0.0
        for n, coefficient in enumerate(PolynomialCoefficients):
            sum += coefficient * np.power(np.log(energy), n)
        sigma = np.exp(sum)
        return sigma

    def rateCoefficientCold(self, crossSection):
        # temperature in eV, sigma in cm^2
        # k=<σv> in m^3/s
        # Assumes a Maxwellian distribution and the neutrals are stationary, projectile energy (as opposed to CoM)
        # I am a bit unsure of using mass vs reducedMass
        mass = crossSection.particle.mass * MP # kg
        sigma = 1e-4 * crossSection.sigma # Convert to m^2
        temperature = self.temperature * QE # Convert to J
        thermalMachNumber = self.MachNumber * np.sqrt(crossSection.particle.charge / 2) # Assumes that Te = Ti

        velocityCOM = np.sqrt(2 * self.energy * QE / (mass)) # Convert to COM, m/s
        thermalVelocity = np.sqrt(2 * temperature / mass)

        if isinstance(self.MachNumber, int) or isinstance(self.MachNumber, float):
            rateCoefficient = np.zeros(len(self.temperature))
            for i, v_th in enumerate(thermalVelocity):
                u = velocityCOM / v_th
                rateCoefficient[i] = v_th / (thermalMachNumber * np.sqrt(np.pi)) * np.trapz(u**2 * sigma * (np.exp(-(thermalMachNumber - u)**2) - np.exp(-(thermalMachNumber + u)**2)), x=u)

        else:
            rateCoefficient = np.zeros((len(self.MachNumber), len(self.temperature)))
            for i, M in enumerate(thermalMachNumber):
                for j, v_th in enumerate(thermalVelocity):
                    u = velocityCOM / v_th
                    rateCoefficient[i, j] = v_th / (M * np.sqrt(np.pi)) * np.trapz(u**2 * (np.exp(-(M - u)**2) - np.exp(-(M + u)**2)) * sigma, x=u)

        rateCoefficient *= 1e6 # m^3/s

        return RateCoefficient(rateCoefficient, self.temperature, crossSection)

    def rateCoefficientHot(self, crossSection):
        # E and T in eV, sigma in cm^2
        # k=<σv> in m^3/s
        # mu, mr is in proton masses
        # Assumes a Maxwellian distribution and the neutrals are hot, COM energy (as opposed to incident)
        # I am a bit unsure of using mass vs reducedMass
        reducedMass = crossSection.reducedMass * MP # kg
        mass = crossSection.particle.mass * MP # kg
        sigma = 1e-4 * crossSection.sigma # Convert to m^2
        temperature = self.temperature * QE # Convert to J
        energy = crossSection.energy * QE # Convert to J

        rateCoefficient = np.array([4 / np.sqrt(2 * np.pi * mass * T) / T * np.trapz(energy * sigma * np.exp(-energy / T), x=energy) for T in temperature])

        rateCoefficient *= 1e6 # cm^3/s

        return RateCoefficient(rateCoefficient, self.temperature, crossSection)
