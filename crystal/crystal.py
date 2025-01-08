import numpy as np
import os
from pathlib import Path
import ase
from typing import Optional, Union, Sequence


class Crystal(ase.Atoms):
    '''
    '''

    def __init__(
            self,
            atoms: ase.Atoms = None,
            pos: np.ndarray = None,
            z_numbers: np.ndarray = None,
            cell: Union[Sequence, np.ndarray] = None,
            occupancy: Union[Sequence, np.ndarray] = None,
        ) -> None:
        '''
        Crystal class for manipulating user structure

        Parameters
        --------
        atoms:  ase.Atoms
            ase.Atoms instance for constructing structure easily.
        
        Returns
        --------
        None
        '''
        self.atoms = atoms if isinstance(atoms, ase.Atoms) else None
        if self.atoms == None and \
            (pos != None & z_numbers != None & cell != None):
            self.pos = np.asarray(pos)
            self.z_numbers = z_numbers
            self.cell = np.eye(cell) if len(cell) == 3 else 
            self.occupancy = occupancy
            self.atoms = ase.Atoms(positions=pos, cell=cell, numbers=z_numbers)
            self.va, self.vb, self.vc,\
                self.alpha, self.beta, self.gamma = self.atoms.get_cell_lengths_and_angles()
            
            self.inv_cell = self.atoms.get_reciprocal_cell()
        else:
            raise RuntimeError('Invalid structure input.')

    @classmethod
    def from_file(cls, path:str):
        '''
        call ase.io.load to read in user structure and return
        an instance of `Crystal`. All supported format according to
        ase.
        '''
        if Path(path).exists():
            return cls(atoms=ase.io.load(path)) 
        else:
            raise ValueError('Invalid path.')

    def calculate_scattering_factor(
            self,
            k_max: float = 2.,
            tol_structure_factor: float = 1e-4,
            return_intensity: bool = False,
    ) -> np.ndarray:
        '''
        calculate scattering intensity
        '''

def xray_scattering_factor(
    s: float,
    a: Union[list[float], np.ndarray],
    b: Union[list[float], np.ndarray],
    c: float
) -> np.ndarray:
    '''
    calculate xray scattering factor

    Args:
        s: float
            scattering vector magnitude \sin(\theta)/\lambda
        a: Union[list[float], np.ndarray]
            Cromer-Mann coefficients \a_i
        b: Union[list[float], np.ndarray]
            Cromer-Mann coefficients \b_i
        c: float
            Cromer-Mann coefficients \c
    
    Returns:
        np.ndarray
            xray scattering factor for an atom
    '''
    # check input
    if not isinstance(a, (list, np.ndarray)):
        raise ValueError('a must be a list or numpy array')
    if not isinstance(b, (list, np.ndarray)):
        raise ValueError('b must be a list or numpy array')
    # change input to numpy array
    a = np.asarray(a)
    b = np.asarray(b)
    # calculate scattering factor
    return np.sum(a * np.exp(-b * s**2)) + c

def electron_scattering_factor(
    s: float,
    a: Union[list[float], np.ndarray],
    b: Union[list[float], np.ndarray],
    c: float,
    Z: int
) -> np.ndarray:
    '''
    calculate electron scattering factor
    
    Args:
        s: float
            scattering vector magnitude \sin(\theta)/\lambda
        a: Union[list[float], np.ndarray]
            Cromer-Mann coefficients \a_i
        b: Union[list[float], np.ndarray]
            Cromer-Mann coefficients \b_i
        c: float
            Cromer-Mann coefficients \c
        Z: int
            atomic number
    
    Returns:
        np.ndarray
            electron scattering factor for an atom
    '''
    # check input
    if not isinstance(a, (list, np.ndarray)):
        raise ValueError('a must be a list or numpy array')
    if not isinstance(b, (list, np.ndarray)):
        raise ValueError('b must be a list or numpy array')
    # change input to numpy array
    a = np.asarray(a)
    b = np.asarray(b)
    # calculate scattering factor
    eta = np.pi * 1e-10 * Z # scaling factor for electron scattering factor
    return xray_scattering_factor(s, a, b, c) / (s**2 + eta**2)

def strucutre_factor(
    atoms: ase.Atoms,
    ) -> np.ndarray:
    
