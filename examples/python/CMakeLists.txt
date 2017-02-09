# Handle installation of the Python module by relying on disttools
install(CODE
  "execute_process(COMMAND python setup.py install -f
  --prefix=${CMAKE_INSTALL_PREFIX}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})")
