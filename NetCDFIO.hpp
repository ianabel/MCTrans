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
		void AddScalarVariable( std::string&& name, std::string&& description, std::string&& units, double value );

	private:
		std::string filename;
		netCDF::NcFile data_file;
};


#endif // NETCDFIO_HPP
