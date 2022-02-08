
#include "Species.hpp"

Species Electron{ .type = Species::Electron, .Charge = -1, .Mass = ElectronMass,  .Name = "Electron" };
Species Proton{ .type = Species::Ion, .Charge = 1, .Mass = ProtonMass, .Name = "Proton"};
Species Deuteron{ .type = Species::Ion, .Charge = 1, .Mass = 1.999*ProtonMass, .Name = "Deuteron" };
Species NeutralHydrogen{ .type = Species::Neutral, .Charge = 0, .Mass = ProtonMass + ElectronMass, .Name = "Neutral Hydrogen"};


