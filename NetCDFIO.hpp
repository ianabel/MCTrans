#ifndef NETCDFIO_HPP
#define NETCDFIO_HPP
#include <netcdf>

/*
 *
 * Class for storing plasma details in a NetCDF file.
 */

class NetCDFIO
{
	public:
		NetCDFIO();
		~NetCDFIO();
		void Open( const std::string &file );
		void Close();
		void AddScalarVariable( std::string   name, std::string description, std::string units, double value );
		void AddTextVariable( std::string name, std::string description, std::string units, std::string text );
		void AddTimeSeries( std::string name, std::string description, std::string units, double initialValue );

		size_t AddTimeSlice( double T );
		void AppendToTimeSeries( std::string const& name, double value, size_t tIndex );

	private:
		std::string filename;
		netCDF::NcFile data_file;
		netCDF::NcDim TimeDim;
		netCDF::NcVar TimeVar;
};


#endif // NETCDFIO_HPP
