SET( ITK_LIBRARIES ${ITK_LIBRARIES} PARENT_SCOPE )

SET(ENABLE_GOOGLE_TEST ${ENABLE_GOOGLE_TEST} PARENT_SCOPE)
IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
