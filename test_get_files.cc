#include "utility.h"

int main( int argc, char **argv ) {
	char **files;
	size_t i;

	if( argc != 2 ) {
		fprintf( stderr, "Usage: %s DIRECTORY\n", argv[0] );
		return EXIT_FAILURE;
	}

	files = get_all_files( argv[1] );
	if( !files ) {
	  fprintf( stderr, "%s: %s: something went wrong\n", argv[0], argv[1] );
	  return EXIT_FAILURE;
	}

	for( i = 0; files[i]; ++i ) {
		if( files[i][0] != '.' ) {
			printf( "%s\n", files[i] );
		}
	}

	for( i = 0; files[i]; ++i ) {
		free( files[i] );
	}
	free( files );

	return EXIT_SUCCESS;
}
