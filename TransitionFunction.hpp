
// Transition Function that is smooth,
// equal to 1.0 for x < L and 0.0 for x > U
// and takes values in (0.0,1.0) for x in (L,U)
double Transition( double x, double L, double U )
{
	if ( L > U ) throw std::invalid_argument( "WHARGL" );
	if ( x <= L ) return 1.0;
	if ( x >= U ) return 0.0;

	double y = ( x - L )/( U - L );
	double arg = 1.0 - 1.0 / ( 1.0 - y*y );
	return ::exp( arg );
}
