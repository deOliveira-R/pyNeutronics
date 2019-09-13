class Isotope:
    def __init__(self, name, symbol, protons, nucleons, mass):
        self.name = name
        self.symbol = symbol
        self.protons = protons
        self.nucleons = nucleons
        self.mass = mass

        # decay_constant = float
        # decay_chain = something...

# class Element:
#     def __init__(self, name, symbol, protons):
#         self.name = name
#         self.symbol = symbol
#         self.protons = protons

H1   =  Isotope('Hydrogen',  'H',  1,   1,   1.00782503223)
O16  =  Isotope(  'Oxygen',  'O',  8,  16,  15.99491461957)
Xe135 = Isotope(   'Xenon', 'Xe', 54, 135, 134.9072278)
Sm149 = Isotope('Samarium', 'Sm', 62, 149, 148.9171921)
U235 =  Isotope( 'Uranium',  'U', 92, 235, 235.0439301)
U238 =  Isotope( 'Uranium',  'U', 92, 238, 238.0507884)
