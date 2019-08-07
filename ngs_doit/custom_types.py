from typing import Optional, List, Tuple, Union, Any, Dict, Callable
from typing_extensions import TypedDict, Literal 

from mypy_extensions import (Arg, DefaultArg, NamedArg,
                             DefaultNamedArg, VarArg, KwArg)
Args = TypedDict('Args', 
        { 'config' : Optional[str],
          'outdir' : Optional[str],
          'ref'    : Optional[str],
          'R1'     : Optional[str],
          'R2'     : Optional[str],
          'sample' : Optional[str] })

PathLike = str
Targets = List[PathLike]
FileDeps = List[PathLike]
Failure = Literal[False]
ActionReturn = Union[Failure, None, Dict[str, Any]]
PyAction = Callable[[List[PathLike], List[PathLike], VarArg(Any)], ActionReturn] 
Actions = List[Union[PyAction, str]]

IdxBamJob = TypedDict('IdxBamJob', 
        { 'targets'  : Targets,
          'actions'  : Actions,
          'file_dep' : FileDeps,
          'task_dep' : List[Literal['index_bam']] })

TaskDepJob = TypedDict('TaskDepJob', 
        { 'targets'  : Targets,
          'actions'  : Actions,
          'file_dep' : FileDeps,
          'task_dep' : List[str] })

Job = TypedDict('Job', 
        { 'targets' :  Targets,
          'file_dep' : FileDeps,
          'actions' :  Actions })

