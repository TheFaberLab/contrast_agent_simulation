# Adding test
option(ENABLE_TESTING "Enable the tests" ${PROJECT_IS_TOP_LEVEL} )

# Note: by default ENABLE_DEVELOPER_MODE is True
# This means that all analysis (sanitizers, static analysis)
# is enabled and all warnings are treated as errors
# if you want to switch this behavior, change TRUE to FALSE
#set(ENABLE_DEVELOPER_MODE
#    TRUE
#    CACHE BOOL "Enable 'developer mode'")
option( ENABLE_DEVELOPER_MODE
        "Enable 'developer mode'"
        ON 
      )
# Change this to false if you want to disable warnings_as_errors in developer mode
option( OPT_WARNINGS_AS_ERRORS_DEVELOPER_DEFAULT
        "Change this to false if you want to disable warnings_as_errors in developer mode"
        ON )



# defaulted_project_options sets recommended defaults and provides user and developer
# modes and full GUI support for choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden
# set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode
# set(<feature_name>_USER_DEFAULT <value>) - set default for user mode
set(ENABLE_CACHE_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(OPT_ENABLE_CLANG_TIDY OFF)