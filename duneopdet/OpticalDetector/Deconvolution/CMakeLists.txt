install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
FILE(GLOB text_files *.txt)
install_source( EXTRAS ${fcl_files} )
install_source( EXTRAS ${text_files} )