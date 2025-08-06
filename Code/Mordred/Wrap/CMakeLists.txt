remove_definitions(-DRDKIT_OSMORDRED_BUILD)
rdkit_python_extension(rdMordred
                       rdMordred.cpp
                       DEST Chem
                       LINK_LIBRARIES Mordred)

# add_pytest(pyMordred ${CMAKE_CURRENT_SOURCE_DIR}/testMordred.py)