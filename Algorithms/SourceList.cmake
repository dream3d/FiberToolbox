
set(${PLUGIN_NAME}_Widgets_HDRS "")
set(${PLUGIN_NAME}_Widgets_SRCS "")

set(${PLUGIN_NAME}_Algorithms_HDRS ${${PLUGIN_NAME}_Algorithms_HDRS}
    ${${PLUGIN_NAME}_SOURCE_DIR}/Algorithms/ComputeGradient.h
)

set(${PLUGIN_NAME}_Algorithms_SRCS ${${PLUGIN_NAME}_Algorithms_SRCS}
    ${${PLUGIN_NAME}_SOURCE_DIR}/Algorithms/ComputeGradient.cpp
)


# Organize the Source files for things like Visual Studio and Xcode
cmp_IDE_SOURCE_PROPERTIES( "${PLUGIN_NAME}/Algorithms" "${${PLUGIN_NAME}_Algorithms_HDRS}" "${${PLUGIN_NAME}_Algorithms_SRCS}" "0")



