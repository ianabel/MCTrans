

#include "Centrifugal.hpp"
#include "FusionYield.hpp"

double PromptAlphaLossFraction( Plasma const& plasma, Configuration const& conf );
double AlphaHeating( Plasma const& plasma, Configuration const& conf );
double AlphaPromptLosses( Plasma const& plasma, Configuration const& conf );
double SlowingDownTime( Plasma const& plasma );

