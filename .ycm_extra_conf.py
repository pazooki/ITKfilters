import os
import ycm_core

#From: https://github.com/oblitum/dotfiles/blob/archlinux/.ycm_extra_conf.py
# '-I/home/phc/repository_local/ITKfilters/src',
# '-I/home/phc/repository_local/ITKfilters/src/Denoise',
# '-I/home/phc/repository_local/ITKfilters/src/Skeleton',
# '-I/home/phc/repository_local/ITKfilters/build-release/gtest-src/include',
# '-I/home/phc/repository_local/ITKfilters/build-release',

flags = [
'-x',
'c++',
'-std=c++11',
'-Wall',
'-Wextra',
'-pedantic',
'-isystem', '/usr/include/c++/6.1.1',
'-isystem', '/usr/include/c++/6.1.1/x86_64-unknown-linux-gnu',
'-isystem', '/usr/include/c++/6.1.1/backward',
'-isystem', '/usr/local/include',
'-isystem', '/usr/include',
'-isystem', '/usr/lib/clang/3.7.1/include',
# '-isystem', '/usr/lib64/gcc/x86_64-unknown-linux-gnu/6.1.1'
'-march=x86-64',
'-mtune=generic',
'-O2',
'-pipe',
'-fstack-protector-strong',
'-frounding-math',
'-msse2',
'-fopenmp',
# '-O3',
'-DNDEBUG',
'-fPIC'
'-DvtkRenderingCore_AUTOINIT=3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL)',
'-isystem', '/home/phc/devtoolset/release/include',
'-isystem', '/usr/include/freetype2',
'-isystem', '/opt/cuda/include/CL',
'-isystem', '/usr/include/GraphicsMagick',
'-isystem', '/usr/include/cairo',
'-isystem', '/usr/include/eigen3',
'-isystem', '/usr/include/qt',
'-isystem', '/usr/include/qt/QtWidgets',
'-isystem', '/usr/include/qt/QtGui',
'-isystem', '/usr/include/qt/QtCore',
'-isystem', '/usr/lib/qt/mkspecs/linux-g++',
'-isystem', '/usr/include/qt/QtOpenGL',
'-isystem', '/usr/include/qt/QtXml',
'-I','./src',
'-I','./src/Denoise',
'-I','./src/Skeleton',
'-I','./src/Quadrature',
'-I','./src/Wavelet',
'-I','./src/Wavelet/IsotropicWaveletFrequencyFunctions',
'-I','./src/Common',
'-I','./build-release/gtest-src/include',
'-I','./build-release',
'-I', '/usr/include/vtk',
'-I', '/home/phc/Software/ITK/install-development/include/ITK-4.10',
# '-I', '/home/phc/Software/ITK/build-development/my-install/include/ITK-4.10',
# '-I/home/phc/devtoolset/release/ITK/include/ITK-4.9',
# '-I/home/phc/devtoolset/release/VTK/include',
]
#Template declaration in .ih files:
#From: https://github.com/Valloric/YouCompleteMe/issues/1938#issuecomment-175411991
# def GetIncludeTemplateDeclaration( filename ):
#     root, ext = os.path.splitext( filename )
#     if ext == '.ih':
#         return [ '-include', root + '.h' ]
#     return []
#
#
# def FlagsForFile( filename, **kwargs ):
#     flags.extend( GetIncludeTemplateDeclaration( filename ) )
#     #Prepend: (no changes)
#     #flags[:0] = GetIncludeTemplateDeclaration( filename )
#
#     return {
#         'flags': flags,
#         'do_cache': True
#     }

# Set this to the absolute path to the folder (NOT the file!) containing the
# compile_commands.json file to use that instead of 'flags'. See here for
# more details: http://clang.llvm.org/docs/JSONCompilationDatabase.html
#
# Most projects will NOT need to set this to anything; you can just change the
# 'flags' list of compilation flags. Notice that YCM itself uses that approach.
# compilation_database_folder = '/home/phc/Software/DGtal/fork-build'
compilation_database_folder = ''

if os.path.exists( compilation_database_folder ):
  database = ycm_core.CompilationDatabase( compilation_database_folder )
else:
  database = None

SOURCE_EXTENSIONS = [ '.cpp', '.cxx', '.cc', '.c', '.m', '.mm' ]


def DirectoryOfThisScript():
  return os.path.dirname( os.path.abspath( __file__ ) )


def MakeRelativePathsInFlagsAbsolute( flags, working_directory ):
  if not working_directory:
    return list( flags )
  new_flags = []
  make_next_absolute = False
  path_flags = [ '-isystem', '-I', '-iquote', '--sysroot=' ]
  for flag in flags:
    new_flag = flag

    if make_next_absolute:
      make_next_absolute = False
      if not flag.startswith( '/' ):
        new_flag = os.path.join( working_directory, flag )

    for path_flag in path_flags:
      if flag == path_flag:
        make_next_absolute = True
        break

      if flag.startswith( path_flag ):
        path = flag[ len( path_flag ): ]
        new_flag = path_flag + os.path.join( working_directory, path )
        break

    if new_flag:
      new_flags.append( new_flag )
  return new_flags

def GetIncludeTemplateDeclaration( filename ):
    root, ext = os.path.splitext( filename )
    if ext == '.hxx':
        return [ '-include', root + '.h' ]
    return []

def FlagsForFile( filename, **kwargs ):
  relative_to = DirectoryOfThisScript()
  final_flags = MakeRelativePathsInFlagsAbsolute( flags, relative_to )
  final_flags.extend(GetIncludeTemplateDeclaration(filename))

  return {
    'flags': final_flags,
    'do_cache': True
  }
