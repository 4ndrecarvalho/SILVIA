# SILVIA project
Modeling the synthesis of SILica materials VIA multiscale computer simulation project. 

This is a silica reative model (__RSi__), based on the MARTINI (Coarse-Grained) Forcefield, that uses a new paradigm for Molecular Dynamics bonding that I developed when I was messing around with Lennard-Jones potentials (LJ). Basically we use dummy atoms to form a intermolecular LJ attractive potential, and use virtual sites that only have the repulsion part of the LJ to impose a well defined geometry for bonding/coordination. This form a small but strong interation pocket that mimics a real bonding interation, and the vibration of the structure will modulate the lability of the bond.

This new reative model applied to CG silica enables the simulation of silica condensation, and formation of silica mesophases in the presence of surfactants in aqueous solutions with excelent agreement with experimental conditions. There are many configurations that you can make, so you can change bead parameters for charge anisotropy and even pH. 

__This was developed for silica but this reative model can be applied to any chemical specie and to any MD package that uses Lennard-Jones potentials (basically all... Amber, LAMMPS, GROMACS, etc...). This model was also tested (initial prototypes) for all atom simulations of Metal Organic Frameworks (MOFs) (link for the video bellow) and proteins within GROMACS.__

https://drive.google.com/file/d/1T85VL0xQN7KZ0D8jxwkJWf8lAdNWDjd-/view?usp=sharing

![Image of RSi](https://lh4.googleusercontent.com/H4cwNgtLouJMLTL3pf9gbUnOkOEbTYogVvDlNna7NLNJrFDNgWFmcrl5TTY9yl7UzhTujVcWuF6YeFlMOouSwSiUX8VM6NtlesSOeXI44HEQS_GfbjGPHQt1s-A5Bgi7lDYTzyCC)


![Image of ENQF21](https://drive.google.com/uc?export=view&id=1PUp9EQEsL0HEHHXFrScnyJlQdio3jnx9)

https://drive.google.com/uc?export=view&id=1PUp9EQEsL0HEHHXFrScnyJlQdio3jnx9


https://drive.google.com/uc?export=view&id=1aLPeQASgzSJtNs9f0Q2Bp7S65RxCkAHr

Molecular dynamics of an aqueous solution of silica (silica condensation).

https://drive.google.com/uc?export=view&id=1DM5H8HeRBwdosoSqQhJWL8vfn6tw5yIi

Molecular dynamics of a silica + CTA+ (cetrilammonium) surfactant in water, with the formation of a MCM41 like mesophase.


https://drive.google.com/uc?export=view&id=1IMysLJ286d_ToRn_yZBKDOB5bsmxh4kw

Details of the final form of the mesophase, with the bigger silica chains highlighted in pink.
 
  
__For more informations read the article:__

[![image](https://user-images.githubusercontent.com/42943782/168496186-989dcd4a-1c8f-4a61-9ba8-2e82a8126973.png)](https://www.nature.com/articles/s41524-022-00722-w)
 
  
CENTRO-01-0145-FEDER-31002 e (PTDC/QUI-QFI/31002/2017).
This work was developed within the scope of project SILVIA, refs. CENTRO-01-0145-FEDER-31002 and PTDC/QUI-QFI/31002/2017 financed by the Portuguese Fundação para a Ciência e a Tecnologia (FCT/MCTES) and co-financed by the European Regional Development Fund (FEDER) under the PT2020 Partnership Agreement.
![image](https://user-images.githubusercontent.com/42943782/168495665-a4bc6a0b-e665-47b6-8415-154f603066ee.png)
![image](https://user-images.githubusercontent.com/42943782/168495672-ca8b6802-035f-4ef4-b7dd-f7fbaaf526e7.png)



