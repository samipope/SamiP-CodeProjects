from __future__ import annotations
import rdkit.Chem
import typing
__all__ = ['PatternFingerprintTautomerTarget', 'TautomerQuery', 'TautomerQueryCanSerialize', 'UnsignedLong_Vect']
class TautomerQuery(Boost.Python.instance):
    """
    The Tautomer Query Class.
      Creates a query that enables structure search accounting for matching of
      Tautomeric forms
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetModifiedAtoms(self) -> UnsignedLong_Vect:
        """
            C++ signature :
                std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> GetModifiedAtoms(RDKit::TautomerQuery {lvalue})
        """
    def GetModifiedBonds(self) -> UnsignedLong_Vect:
        """
            C++ signature :
                std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> GetModifiedBonds(RDKit::TautomerQuery {lvalue})
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatch(RDKit::TautomerQuery,RDKit::ROMol [,bool=False [,bool=False]])
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatch(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatches(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatches(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])
        """
    @typing.overload
    def GetSubstructMatches(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatches(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatchesWithTautomers(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])
        """
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            C++ signature :
                _object* GetSubstructMatchesWithTautomers(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    def GetTautomers(self) -> typing.Any:
        """
            C++ signature :
                _object* GetTautomers(RDKit::TautomerQuery)
        """
    def GetTemplateMolecule(self) -> rdkit.Chem.Mol:
        """
            C++ signature :
                RDKit::ROMol GetTemplateMolecule(RDKit::TautomerQuery {lvalue})
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
            C++ signature :
                bool IsSubstructOf(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False]]])
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, params: SubstructMatchParameters) -> bool:
        """
            C++ signature :
                bool IsSubstructOf(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    def PatternFingerprintTemplate(self, fingerprintSize: int = 2048) -> ExplicitBitVect:
        """
            C++ signature :
                ExplicitBitVect* PatternFingerprintTemplate(RDKit::TautomerQuery {lvalue} [,unsigned int=2048])
        """
    def ToBinary(self) -> typing.Any:
        """
            C++ signature :
                boost::python::api::object ToBinary(RDKit::TautomerQuery)
        """
    def __getinitargs__(self) -> tuple:
        """
            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::TautomerQuery)
        """
    def __getstate__(self) -> tuple:
        """
            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @typing.overload
    def __init__(self, arg1: Mol) -> typing.Any:
        """
            C++ signature :
                void* __init__(boost::python::api::object,RDKit::ROMol)
        """
    @typing.overload
    def __init__(self, arg1: Mol, arg2: str) -> typing.Any:
        """
            C++ signature :
                void* __init__(boost::python::api::object,RDKit::ROMol,std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>)
        """
    @typing.overload
    def __init__(self, pickle: str) -> None:
        """
            C++ signature :
                void __init__(_object*,std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>)
        """
    def __setstate__(self, data: tuple) -> None:
        """
            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
class UnsignedLong_Vect(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
            C++ signature :
                bool __contains__(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue},_object*)
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
            C++ signature :
                void __delitem__(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue},_object*)
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::__1::vector<unsigned long, std::__1::allocator<unsigned long>>&>,_object*)
        """
    def __init__(self) -> None:
        """
            C++ signature :
                void __init__(_object*)
        """
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__wrap_iter<unsigned long*>> __iter__(boost::python::back_reference<std::__1::vector<unsigned long, std::__1::allocator<unsigned long>>&>)
        """
    def __len__(self) -> int:
        """
            C++ signature :
                unsigned long __len__(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue})
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
            C++ signature :
                void __setitem__(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue},_object*,_object*)
        """
    def append(self, item: typing.Any) -> None:
        """
            C++ signature :
                void append(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue},boost::python::api::object)
        """
    def extend(self, other: typing.Any) -> None:
        """
            C++ signature :
                void extend(std::__1::vector<unsigned long, std::__1::allocator<unsigned long>> {lvalue},boost::python::api::object)
        """
def PatternFingerprintTautomerTarget(target: Mol, fingerprintSize: int = 2048) -> ExplicitBitVect:
    """
        C++ signature :
            ExplicitBitVect* PatternFingerprintTautomerTarget(RDKit::ROMol [,unsigned int=2048])
    """
def TautomerQueryCanSerialize() -> bool:
    """
        Returns True if the TautomerQuery is serializable (requires that the RDKit was built with boost::serialization)
    
        C++ signature :
            bool TautomerQueryCanSerialize()
    """
