# - Try to find ECS    
# Once done this will define
#
#  ECS_EXE          - ECS Path       
#

set(ECS_FOUND FALSE)

# Check for ecs_join
find_program(ECS_EXE
  NAMES ecs_join  
  HINTS  $ENV{ECSPATH}
  NO_DEFAULT_PATH
  DOC "The prepro saturne executable"
  )
find_program(ECS_EXE
  NAMES ecs_join   
  DOC "The prepro saturne executable"
  )

# Standard package handling
find_package_handle_standard_args(ECS
                                  "ECS could not be found. Be sure to set ECSPATH environment variable."
                                  ECS_EXE)

mark_as_advanced(ECS_EXE)


