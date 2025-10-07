# ^Modeling the Cleavage Efficiency of BoNT/A Light Chain^

^1. Objective^
![^cap^](https://static.igem.wiki/teams/5924/drylab/tail.webp)
    When botulinum toxin cleaves SNAP25, it leaves behind a 13-peptide tail, and since the alpha factor is also a 13-peptide in size, this raises a question: would the remaining tail affect the binding of the alpha factor to Ste2? If so, how can the SNAP25 sequence be modified to ensure that, while maintaining sufficient cleavage activity by botulinum toxin, the tailed alpha factor can still bind to Ste2 and trigger downstream responses?

    To study the differences in cleavage efficiency of the BoNT/A light chain on SNAP-25 mutants (specifically, deletions in the 199-202 region) and to combine this analysis with the binding affinity of the resulting cleavage products to the Ste2 receptor. The ultimate goal is to screen for the optimal mutant sequence.

^2. Major Challenges^

     ^Scarcity of Experimental Enzyme Kinetics Data^

        Initially, the plan was to quantify the impact of different amino acid mutations (in the 199-202 region) on cleavage efficiency by using corresponding enzyme kinetic parameters (kcat/Km). However, preliminary research revealed that such experimental data is extremely scarce, making this approach unfeasible.

     ^Inadequate Adaptability of Deep Learning Models for Protein Substrates^

        Subsequently, data-driven deep learning models were considered for predicting enzyme kinetic parameters. Benchmark tools like the CataPro model (Wang et al., 2025) are powerful but do not support reactions where the substrate is a protein. Fine-tuning a new model from scratch would be too time-consuming and would be hampered by a lack of sufficient training data.

![](https://static.igem.wiki/teams/5924/drylab/catapro.webp)

        Figure 1: CataPro for predicting enzyme kinetic parameters of enzymatic reactions.

     ^Need to Quantify Both Cleavage Efficiency and Downstream Binding Effects^

        A comprehensive solution requires a method that can simultaneously quantify both the efficiency of the enzymatic cleavage and the binding affinity of the subsequent products to their downstream targets.

^3. Calculation of Binding Efficiency^

	 ^Structure Prediction and Preprocessing ^

        a.AlphaFold2 Prediction:

        Input:The amino acid sequences of the constructed full-length protein and four deletion mutants (in the 199-202 region).

        Output: Five 3D protein structures (assessing confident regions with pLDDT > 80).

        b.Structure Optimization:

        Loop regions were optimized using Modeller.

        Energy minimization was performed using AmberTools.

     ^Molecular Binding Efficiency Analysis^

        a.Inputs:

        Crystal structure of the BoNT/A light chain.

        Predicted structures of the SNAP-25 mutants.

        b.Tool:

        ClusPro 2.0, which performs rigid-body docking based on the Fast Fourier Transform (FFT) algorithm.

        c.Docking Process:

        ClusPro 2.0 first generates 70,000 different rotations of the ligand (e.g., SNAP-25 mutant). For each rotation, it performs translations along the x, y, and z axes relative to the receptor (e.g., BoNT/A light chain crystal structure), selecting the best-scoring translational combination for that rotation. From the 70,000 initial combinations, the 1,000 with the lowest scores are selected. These 1,000 conformations are then subjected to a greedy clustering analysis with a 9 Å C-alpha RMSD radius to identify the most likely binding poses.

        d.Docking Result:

        The output from ClusPro 2.0 typically consists of 1,000 low-energy results, which are clustered to find the native binding site with the widest free energy funnels. Results are presented in various formats, including the docked complex structures and scores for each model.
![Docking results of BoNT/A with SNAP-25 variants.](https://static.igem.wiki/teams/5924/drylab/snap-25.webp)

        Figure 2: Docking results of BoNT/A with SNAP-25 variants. 1: Wild-type protein. 2: Deletion of residue 202. 3: Deletion of 201-202. 4: Deletion of 200-202. 5: Deletion of 199-202. A lower score indicates higher potential cleavage efficiency.

![Docking results of cleavage products with Ste2 receptor.](https://static.igem.wiki/teams/5924/drylab/ste2.webp)

        Figure 3: Labels correspond to the variants in Figure 2. These are the docking results of the different cleavage products with the Ste2 receptor protein. A lower score indicates stronger binding affinity. 6 is the control group, measuring the binding affinity between alpha-factor and the Ste2 receptor.

    ^Results and Validation (Revised Method)^

        To obtain a more robust and physically meaningful metric for binding affinity, the initial approach of using a simple weighted average has been replaced with a method based on statistical mechanics. This approach correctly accounts for the contributions of all significant binding modes (clusters) based on their energy and population size.

    ^Boltzmann-Weighted Score (BWS)^

        Instead of a simple composite score, a Boltzmann-Weighted Score (BWS) is calculated for each construct. This score aggregates the contributions of all docked clusters, giving more weight to lower-energy and more populated clusters, which is physically more representative of the overall binding propensity.  

        The formula is:

        \[text{BWS} = \sum_{i} \left( N_i \times e^{-E_i / kT} \right)\] 

        Parameter Definitions:

        *i*: An index for each cluster generated by ClusPro for a given construct.

        *N_i*: The number of members (population) in cluster i. This reflects the size of the conformational basin for that binding mode.

        *E_i*: The representative energy score of cluster i. It is recommended to use the cluster's Center score, as it is a better representation of the cluster's average energy and is less sensitive to single outliers than the Lowest_Energy score.


        *kT*: The thermal energy term. For ranking purposes with docking scores, kT can be effectively set to 1. This simplifies the calculation while preserving the crucial exponential relationship between energy and probability, ensuring that lower-energy states are weighted more heavily.

        ^Relative Binding Efficiency (RBE)^

        Using the BWS, a direct and meaningful Relative Binding Efficiency (RBE) can be calculated. This ratio reflects how the overall binding propensity of a construct compares to a reference.   

        The formula is:   

        \[text{RBE}_{\text{construct}} = \left( \frac{\text{BWS}_{\text{construct}}}{\text{BWS}_{\text{reference}}} \right) \times 100\%\]

        Parameter Definitions:

        *BWS_construct*: The Boltzmann-Weighted Score calculated for the construct of interest.

        *BWS_reference*: The BWS of the baseline construct (e.g., the Control α-factor or the wild-type SNAP25 partial). This sets the 100% reference point for comparison.

        This revised methodology provides a more accurate and scientifically defensible ranking of the constructs' binding affinities.


^4. Modeling the Dynamic Cleavage and Product Release Process^

    Enzymatic cleavage is inherently a dynamic event. Standard static docking is insufficient as it only captures the initial binding pose and cannot model the subsequent transition state, bond scission, and product release. To address this, a more sophisticated multi-step simulation approach is necessary to effectively model the process and estimate cleavage efficiency.

    ^Defining the Post-Cleavage State^

        The first step is to accurately represent the system immediately after the peptide bond is broken.

    ^Key Operations^

        Substrate Cleavage Site Definition:

        Based on the experimentally confirmed BoNT/A cleavage site (Q197-R198 of SNAP-25), the substrate structure is computationally modified to generate two distinct product fragments (an N-terminal and a C-terminal peptide). Instead of an abrupt covalent break, a more physically realistic approach is used, employing the nmropt function in AMBER to gradually weaken the force constant of the target peptide bond, simulating a bond-softening and dissociation pathway.

    ^Modeling the Enzyme-Product Complex with Flexible Docking^

        The interaction between the enzyme and the two newly formed product fragments is then modeled to understand the product release energetics.

    ^Key Steps^

        a.Docking Software:

         HADDOCK 2.4 is an excellent choice, selected for its strength in handling multi-body systems and its ability to incorporate data-driven restraints.

        b.Flexible Region Definition: 
        
        To account for the conformational rearrangements that occur upon cleavage and during product release, crucial regions must be defined as flexible. These include:

            The catalytic active site of the BoNT/A light chain (LC).

            The residues surrounding the newly formed C- and N-termini of the two SNAP-25 product fragments.

        c.Application of Restraints:
        
        To guide the simulation towards a realistic dissociated state, distance restraints can be applied.

        An ambiguous restraint ensuring the catalytic residues of BoNT/A remain proximal to the cleavage site.

        A negative restraint (e.g., distance > 8 Å) between the new C-terminus (of the first fragment) and N-terminus (of the second fragment) can be used to model their separation and prevent unrealistic re-ligation.

    ^Conformational Sampling and Selection:^
    
        A large ensemble of complex conformations (e.g., 1,000) is generated.

        These are scored and ranked (e.g., using the HADDOCK score), and the top-ranking models are selected.

        Cluster analysis (e.g., based on RMSD with a threshold of 5.0 Å) is performed to identify the most representative and populated conformational states of the enzyme-product complex.  

    ^Results and Interpretation (Post-Simulation Analysis)^

        The energy landscapes derived from these more advanced simulations, often represented by calculating the activation energy barriers for product release (as conceptualized in Figure 6), provide a more direct and physically relevant proxy for cleavage efficiency than simple binding scores. By comparing these energy barriers across different mutants, a robust ranking of their relative cleavage efficiencies can be established. Based on a comprehensive analysis, the 199-202Δ mutant can be identified as the most suitable candidate among the tested proteins.

![](https://static.igem.wiki/teams/5924/drylab/cleavage.webp)

    Figure 6: Activation energy calculation results, representing a relative measure of cleavage efficiency.

 ^5. Conclusion^

    This study aimed to identify an optimal SNAP-25 mutant by evaluating both its cleavage efficiency by the BoNT/A light chain and the subsequent binding affinity of its product to the Ste2 receptor. Faced with a lack of experimental kinetic data, we developed a comprehensive in silico workflow combining structure prediction, multi-stage molecular docking, and a physically-grounded scoring methodology based on Boltzmann statistics.

    Our integrated analysis consistently identifies the SNAP-25 partialΔ199-202 mutant as the most promising candidate. This conclusion is supported by two key findings:

    ^Enhanced Cleavage Efficiency^: 
    
    Dynamic modeling of the product release process indicated that the Δ199-202 mutant possesses the lowest activation energy barrier, suggesting a faster and more efficient catalytic turnover.

    ^Favorable Downstream Binding^:
    
    The cleavage product of this mutant maintained a strong binding affinity for the Ste2 receptor, ensuring that its subsequent biological function is not compromised.

    The findings provide a strong, data-driven rationale for prioritizing the Δ199-202 mutant for experimental validation. Furthermore, the computational pipeline established here serves as a robust framework for similar enzyme engineering and substrate optimization problems, particularly in data-scarce scenarios. Future work should focus on experimentally verifying these predictions through in vitro assays.

^6. References^

[1]Fei Wang, Honglin Wan, Jianping Hu, and Shuang Chang. Molecular dynamics simulations of wild type and mutants of botulinum neurotoxin a complexed with synaptic vesicle protein 2c. Molecular BioSystems, 11(1):223–231, 2015. doi: 10.1039/c4mb00383g.

[2]Wang, Z., Xie, D., Wu, D., Luo, X., Wang, S., Li, Y., Yang, Y., Li, W., & Zheng, L. (2025). Robust enzyme discovery and engineering with deep learning using CataPro. Nature communications, 16(1), 2736. https://doi.org/10.1038/s41467-025-58038-4


# Molecular dynamics simulation between Ste 2 and the fusion protein

Our biosensor relies on the binding of the α-factor released after toxin cleavage to the Ste2 receptor. A key scientific question is: If the fusion protein(contains Pir1p,Snap25 and α-factor) detaches from the cell wall at a certain probability  ,could the α-factor on the free fusion protein bind to Ste2 , leading to false positive signals?

Using molecular dynamics simulations not only helps us intuitively observe the behavior of the two components, but also provides reliable quantitative data to verify the reliability of the biosensor. Therefore, we conducted molecular dynamics simulations and performed a detailed analysis of the output data. Through simulations exceeding 100 nanoseconds, we confirmed that the steric hindrance effect generated by the linker effectively isolated the α-factor from the Ste2 receptor, eliminating the risk of false positives during the simulation period and providing crucial theoretical support for subsequent wet experiments.

Additionally, our team has provided a quick guide to molecular dynamics simulation, making this tool more accessible and user-friendly.

^1. Preparation for operation: Software environment and hardware equipment^

We chose to use Gromacs as the simulation tool and Pymol as the visualization tool. All molecular dynamics simulations were performed on a high-performance workstation running Ubuntu 20.04.6 LTS, equipped with an Intel Core i9-14900KF CPU (24 cores/32 threads) and an NVIDIA GeForce RTX 4070 GPU (12GB VRAM). This configuration provided substantial parallel processing capabilities through 32 available CPU threads and GPU acceleration via CUDA 12.8, enabling efficient execution of our 100 ns production simulations using GROMACS with full hardware optimization.

^2. System Setup: Molecular Model Construction and Force Field Selection^

The basis of molecular dynamics simulations is the molecular data and force field configuration of the research subject. 
In terms of constructing the research subject, we obtained three key structures constructed by the experimental group: the Ste2 receptor (receptor.pdb), the fusion protein (containing the SNAP-25 linker α factor ,snap25 and  Pir1p (fusion.pdb)), and the Ste2 receptor lacking the α factor (receptor_ref.pdb) as a control.

![From left to right are the fusion protein, the Ste2 receptor with α factor, and the Ste2 receptor without α factor.](https://static.igem.wiki/teams/5924/drylab/molecular.webp)

For the force field and other configurations, we chose the CHARMM36 force field, and for the solution box, we selected the SPCE water box.

```bash
align fusion and chain A, receptor_ref  
select alpha_linker , fusion and resi 407-414
save initial_complex.pdb, fusion
gmx pdb2gmx -f initial_complex.pdb -o processed.gro -p topol.top -water spce -ff charmm36 -ignh  
```

Our goal is to verify whether, in a solution environment, the normal molecular motions of the α linked to the N- and C-termini of the fusion protein might cause unintended binding with STE2, thereby producing false positive signals and interfering with toxin detection. Therefore, based on three key structures, we used PyMOL to spatially align the N-terminus of the fusion protein with the extracellular domain of the receptor, simulating its display state on the yeast cell wall. The generated initial_complex.pdb is the initial structure file:

^3. Build a balanced simulation system: Energy Minimization and NVT & NPT Equilibration^

Energy minimization is the primary step in molecular dynamics simulations, with its core purpose being to eliminate atomic overlaps and geometric conflicts in the initial model. By optimizing atomic positions, it relaxes the system to a lower-energy, stable state, thereby removing unreasonable internal stresses and providing a physically reasonable starting point for subsequent equilibration and production simulations, ensuring numerical stability and reliability of the results.

```bash
# Define a simulation box
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
# Solvation
gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top
# Add particles to physiological concentration
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
# Using the steepest descent method, the energy is converged to a minimum.
gmx grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr 
gmx mdrun -v -deffnm em -nt 1
```

After completing the energy minimization step, the potential energy decreased from the initial -4,447,399.50 kJ/mol to the final -13,915,840.00 kJ/mol, a drop of 212.9%.

![Energy Mnimization](https://static.igem.wiki/teams/5924/drylab/energy-minimization-validation.webp)

After completing energy minimization, allow the system to relax under the set temperature and pressure to reach a stable dynamic equilibrium that corresponds to real experimental conditions (such as specific temperature, pressure, and density), preparing the system for production simulation. This work is divided into two consecutive stages: NVT equilibration (constant number of particles, volume, and temperature) and NPT equilibration (constant number of particles, pressure, and temperature). During this process, positional restraints are applied to the protein atoms.

a. NVT equilibrium, also known as isothermal equilibrium, aims to assign the correct temperature (kinetic energy distribution) to the system. Starting from the 'cold' system (0 K) obtained from energy minimization, the temperature of the system is gradually increased and stabilized at our target temperature of 310 K (37°C, physiological temperature) within 100 ps using a thermostat (we use V-rescale). During this phase, the system volume is fixed, which means the system density will adjust itself as the temperature rises.

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v -nt 2
```

After the NVT equilibration, we obtained a system with the correct temperature but the density may not yet be at the target value.

b. NPT equilibration, also known as isothermal-isobaric equilibration, aims to adjust the system density to the correct value by controlling the pressure. Based on the NVT stage, we maintain the temperature at 310 K. At the same time, we use the Parrinello-Rahman method to enable pressure coupling, keeping the system pressure stable at 1 bar (standard atmospheric pressure). The simulation box can then fluctuate in volume. The pressure coupler will adjust the box size to bring the internal pressure of the system to the target value, allowing the density of water molecules and ions to reach approximately that of real water, ~1000 kg/m³. In this stage, we do not regenerate velocities, but continue using the velocities from the NVT stage, allowing a smooth transition from NVT to NPT.

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v -nt 2
```

![NVT & NPT](https://static.igem.wiki/teams/5924/drylab/nvt-npt-equilibration-validation.webp)

After the NPT equilibration, we obtained a system with correct temperature, pressure, and density, and the solvent molecules had formed a stable environment around the protein.

^4.Production simulation and data analysis^

Based on our molecular dynamics work, we performed a 100-ns production simulation to study the dynamic behavior of the fusion protein and its interaction with the Ste2 receptor. The simulation was carried out using GROMACS under the CHARMM36 force field, with a 2-fs time step under NPT ensemble conditions at 310 K and 1.0 bar. The simulation aimed to observe whether structural fluctuations in the uncut fusion protein would lead to unintended binding—a key question for assessing potential false positives in our biosensor system. The detailed parameter table of our simulation is as follows.

| Parameter Category       | Parameter Name         | Value/Setting       | Description                                                                 |
|---------|----------------|---------------------|-----------------------------------------------------------------------------|
| Simulation Length    | nsteps                 | 50,000,000          | Total steps, corresponding to 100 ns (step size 0.002 ps)                   |
| Time Step            | dt                     | 0.002               | 2 fs                                                                        |
| Coordinate Output    | nstxout                | 0                   | Do not output detailed coordinate trajectories (save space)                 |
| Velocity Output      | nstvout                | 0                   | Do not output velocity trajectories                                         |
| Force Output         | nstfout                | 0                   | Do not output force trajectories                                            |
| Compressed Coord Output | nstxout-compressed  | 5000                | Output compressed coordinates every 10 ps                                   |
| Energy Output        | nstenergy              | 5000                | Output energy every 10 ps                                                   |
| Log Output           | nstlog                 | 5000                | Output log every 10 ps                                                      |
| Bond Constraint Algorithm | constraint_algorithm | lincs             | Use LINCS algorithm for bond constraints                                    |
| Constraint Type      | constraints            | h-bonds             | Constrain all hydrogen-related bonds                                        |
| Neighbor Searching   | cutoff-scheme          | Verlet              | Use Verlet list update scheme                                               |
| Neighbor Search Type | ns_type                | grid                | Use grid search                                                             |
| Neighbor List Update | nstlist                | 20                  | Update neighbor list every 20 steps                                         |
| Electrostatic Cutoff | rcoulomb               | 1.0                 | Short-range electrostatic cutoff radius (nm)                                |
| Van der Waals Cutoff | rvdw                   | 1.0                 | Van der Waals interaction cutoff radius (nm)                                |
| Electrostatic Method | coulombtype            | PME                 | Use Particle Mesh Ewald for long-range electrostatics                       |
| Temperature Coupling | tcoupl                 | V-rescale           | Use V-rescale temperature coupling method                                   |
| Temperature Groups   | tc-grps                | Protein Non-Protein | Couple temperature separately for protein and non-protein groups            |
| Reference Temperature| ref_t                  | 310 310             | Reference temperature 310 K for both groups                                 |
| Pressure Coupling    | pcoupl                 | Parrinello-Rahman   | Use Parrinello-Rahman pressure coupling method                              |
| Pressure Coupling Type | pcoupltype           | isotropic           | Isotropic pressure coupling                                                 |
| Reference Pressure   | ref_p                  | 1.0                 | Reference pressure 1.0 bar                                                  |
| Isothermal Compressibility | compressibility | 4.5e-5           | Isothermal compressibility of water (bar⁻¹)                                 |
| Continuation         | continuation           | yes                 | Continue simulation from equilibrium phase                                  |
| Velocity Generation  | gen_vel                | no                  | Use velocities generated from previous equilibrium phase                    |
| Periodic Boundaries  | pbc                    | xyz                 | Use periodic boundary conditions in x, y, z directions                      |
| Position Restraints  | define                 | (commented out)     | No position restraints on protein in production run                         |

This simulation used GPU acceleration with assistance from a 16-core CPU, and the runtime was 46 hours.

```bash
# Prepare for production simulation input
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Start production simulation
nohup /usr/local/gromacs_gpu/bin/gmx mdrun -s md_100ns.tpr -deffnm md_100ns -v \
    -ntmpi 1 -ntomp 1 \
    -nb gpu -bonded gpu -pme gpu \
    -gpu_id 0 \
    > md_100ns.log 2>&1 &

# Cycle calibration and stacking eliminate cyclical boundary effects.
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_noPBC.xtc -pbc mol -center
```

We conducted the following analysis on the obtained 100 ns simulation results:

a. We directly visualize the trajectories obtained from the production simulations, and the demonstration in Pymol shows that the $alpha$ factor wanders around the ste2 receptor under the traction of the linker but never approaches the binding pocket.

![visualization](https://static.igem.wiki/teams/5924/drylab/movie-final.webp)

(the fusion protein and the ste2 were painted cyan)


b. Turning radius analysis

Radius of gyration analysis can confirm the conformational stability of fusion proteins from the perspective of structural dynamics, providing key molecular-level evidence for the low false-positive risk of biosensors.

```bash
gmx gyrate -s md_5ns.tpr -f md_5ns.xtc -o gyrate.xvg
```

Based on a systematic analysis of a 100 ns molecular dynamics simulation trajectory, the calculation of the radius of gyration revealed that the fusion protein maintains a stable spatial conformation in solution, with an average size of 3.86 nm and relatively small fluctuations (approximately 4.6% relative fluctuation). This indicates that the protein structure is well-packed, without significant folding or extension, further supporting from a structural dynamics perspective that it can effectively avoid unintended interactions with the Ste2 receptor in its intact state, providing key theoretical guidance for designing low false-positive biosensors.

![Turing radius analysis](https://static.igem.wiki/teams/5924/drylab/comprehensive-radius-of-gyration-analysis.webp)

From Table A, it can be seen that the average radius of gyration is 3.856 nm, with a fluctuation range of 3.482 - 4.140 nm (±0.334 nm), a standard deviation of 0.176 nm, and a relative fluctuation of 4.6%. A relative fluctuation of 4.6% indicates that the conformational changes of the fusion protein are within an acceptable range, and the protein maintains a relatively stable tertiary structure. The size of the fusion protein, 3.86 nm, is suitable for display on the cell surface, being neither too compact (which could affect function) nor excessively extended (which could lead to degradation).

As seen in Table B, the radius of gyration in the X direction is 2.355 nm (most compact), in the Y direction is 3.590 nm (most extended), and in the Z direction is 3.357 nm (moderate), with an anisotropy of 1.235 nm. The Y-direction is 52% longer than the X-direction, indicating that the protein is ellipsoidal in shape. This shape is conducive to the exertion of steric effects. Anisotropic shapes can more effectively isolate the α factor, form asymmetric structures, and reduce the probability of accidental contact with the Ste2 receptor.


c. RMSD
We calculated the RMSD of the protein backbone atoms to assess whether the system had reached a stable state during the production simulation.

```bash
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
```

RMSD analysis indicates that the fusion protein remains structurally stable and convergent during the 100 ns simulation, providing key dynamic evidence to confirm that its spatial conformation does not lead to false-positive binding.

![RMSD](https://static.igem.wiki/teams/5924/drylab/rmsd-plot-scaled.webp)

The average RMSD was 0.65 ± 0.25 nm, with a fluctuation range of 0.40 - 1.35 nm. The system reached a stable equilibrium after approximately 40 ns. The RMSD values remained below 1.0 nm, indicating that the fusion protein maintained an intact folded structure. There were no drastic conformational changes, supporting the structural reliability of the sensor.


Through 100-nanosecond all-atom molecular dynamics simulations, we validated at the atomic level that the fusion protein can effectively maintain its structural integrity in the uncleaved state, with the alpha factor and Ste2 receptor maintaining a stable spatial separation, and the theoretical risk of false positives is extremely low. Although our molecular dynamics simulations provide valuable theoretical insights for sensor design, it is important to acknowledge a significant limitation: the simulations did not include the cellular membrane environment of the Ste2 receptor. On the surface of a real yeast cell, membrane fluidity, local concentration effects, and spatial constraints may alter the interaction dynamics between the fusion protein and the receptor. Therefore, we strongly recommend that the experimental team validate the sensor's specificity in the real membrane environment after constructing the yeast strains. 


^5. Reference^

[1] Madej T, Lanczycki CJ, Zhang D, Thiessen PA, Geer RC, Marchler-Bauer A, Bryant SH. " MMDB and VAST+: tracking structural similarities between macromolecular complexes. Nucleic Acids Res. 2014 Jan; 42(Database issue):D297-303

[2]Abramson, J., Adler, J., Dunger, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). https://doi.org/10.1038/s41586-024-07487-w

# ODE for simulation of the fusion protein shedding

^1. Abstract^

We established a system of ordinary differential equations to quantitatively simulate and analyze the expression of the fusion protein Pir1P-SNAP25-α factor in Substrate-acting Strain, its loss during lyophilization, and its dynamic changes upon rehydration. The results indicated that yeast cells rapidly returned to steady-state expression after rehydration, with the number of cell wall–anchored fusion proteins remaining relatively stable, while the proportion of spontaneously detached proteins consistently low. Therefore, we conclude that the spontaneous detachment of the fusion protein has a negligible effect on the accuracy of botulinum toxin detection, ensuring the reliability of the assay system.

^2. Introduction^

In our project, the fusion protein gene was stably integrated into the the genome of the substrate-acting strain to enable continuous expression. The substrate-acting strain is further intended to be preserved by lyophilization, formulated as ready-to-use product, and rehydrated before detection. Based on this workflow, gene expression modeling of the substrate-acting strain was divided into three stages: cultivation, lyophilization, and rehydration.

^3. Model^

 a. Cultivation Stage

To produce our product,the substrate-acting strain with the integrated fusion protein gene was continuously cultured for an extended period before lyophilization. During this cultivation stage, the expression of the fusion protein proceeds sequentially through transcription, translation, protein processing, and secretion and anchoring (as illustrated in Fig. 1) [1,2].

![Fig. 1 The process of fusion protein Pir1P-SNAP25-α factor expression and anchoring](https://static.igem.wiki/teams/5924/drylab/gene-expression-flowchart.webp)

1) Transcription: 

\[
\frac{dM}{dt} = k_{\text{tx}} - (k_{M,\text{deg}} + k_{\text{dil}}) \cdot M
\]

2) Translation (mRNA is translated into nascent, immature polypeptide chains):

\[
\frac{dP_{\text{immature}}}{dt} = k_{\text{tl}} \cdot M - (k_{\text{mature}} + k_{\text{Pc,deg1}} + k_{\text{dil}}) \cdot P_{\text{immature}}
\]

3) Protein Processing (immature peptides mature through folding and post-translational modifications):

\[
\frac{dP_{\text{mature}}}{dt} = k_{\text{mature}} \cdot P_{\text{immature}} - (k_{\text{wall}} + k_{\text{Pc,deg2}} + k_{\text{dil}}) \cdot P_{\text{mature}}
\]

4) Secretion & Anchoring (mature fusion proteins are released into the periplasmic space of yeast and covalently anchored to the yeast cell wall):

\[
\frac{dP_{\text{wall}}}{dt} = k_{\text{wall}} \cdot P_{\text{mature}} - (k_{\text{free}} + k_{\text{dil}}) \cdot P_{\text{wall}}
\]

5) Spontaneous Release (A small proportion of fusion proteins are spontaneously released from the cell wall):

\[
\frac{dP_{\text{free}}}{dt} = k_{\text{free}} \cdot P_{\text{wall}}
\]

^4. Parameter List^

|Parameter |Meaning |
|:---: |:---|
|$M$ |Copy number of mRNA corresponding to the fusion protein |
|$P_{\text{immature}}$ |Number of nascent polypeptides translated from mRNA |
|$P_{\text{mature}}$ |Number of mature fusion proteins |
|$P_{\text{wall}}$ |Number of fusion proteins anchored to the yeast cell wall |
|$P_{\text{free}}$ |Number of fusion proteins released from the cell wall |

|Parameter |Meaning |Value |References |
|:---: |:--- |:---: |:---: |
|$k_{\text{tx}}$ |Translation rate of the fusion protein |0.25/min |estimated, [10] |
|$k_{M,\text{deg}}$ |Rate of RNA degradation |0.0347/min |[11] |
|$k_{\text{dil}}$ |Dilution rate of yeast growth |0.00167/min |[12] |
|$k_{\text{tl}}$ |Rate of translation |0.96/min |[13] |
|$k_{\text{mature}}$ |Rate of protein maturation |0.04621/min |[14] |
|$k_{\text{Pc,deg1}}$ |Rate of peptide degradation |0.0693/min |[15] |
|$k_{\text{Pc,deg2}}$ |Rate of protein degradation |0.01612/min |[16] |
|$k_{\text{wall}}$ |Rate of mature protein translocation and anchoring to the cell wall |0.03/min |[17] |
|$k_{\text{free}}$ |Rate of spontaneous release of the fusion protein |$4.81\times10^{-5}$/min |[18] |


The simulation process throughout the cultivation stage is shown in Fig. 2.

![Fig. 2 Cultivation Stage](https://static.igem.wiki/teams/5924/drylab/cultivation-stage.webp)

Because the cultivation period is sufficiently long, the gene expression process reaches steady state prior to lyophilization. $M$,$P_{\text{immature}}$,$P_{\text{mature}}$,$P_{\text{wall}}$,$P_{\text{free}}$,etc.,are constants,and therefore we define:

1) Transcription
\[
    \frac{dM}{dt} = 0
\]
2) Translation(mRNA translation into nascent, immature polypeptide chains):
\[
    \frac{dP_{\text{immature}}}{dt} = 0
\]
3) Protein Processing(Maturation of immature peptides through folding and post-translational modifications):
\[
    \frac{dP_{\text{mature}}}{dt} = 0
\]
4) Secretion & Anchoring(Release of the mature fusion protein into the periplasmic space of yeast, followed by covalent anchoring to the yeast cell wall):
\[
    \frac{dP_{\text{wall}}}{dt} = 0
\]
At steady state, the molecular abundances were calculated as:
\[
    M = 6.9 molecules/cell
\]
\[
    P_{\text{immature}} = 56.3 molecules/cell
\]
\[
    P_{\text{mature}} = 54.5 molecules/cell
\]
\[
    P_{\text{wall}} = 950.8 molecules/cell
\]

In the cultivation stage, we assumed that once secreted into the periplasmic space of yeast, PIR-fusion proteins were immediately anchored to the cell wall, thereby neglecting any time delay in anchoring. This assumption is justified since PIR proteins can rapidly undergo covalent attachment to the yeast cell wall under the action of specific enzymes [3], and no reliable experimental data for anchoring kinetics are currently available.

 b. Lyophilization Stage:

During lyophilization, multiple biomolecules within and outside yeast cells inevitably suffer structural and functional damage. This may result from low-temperature stress, ice crystal formation, conformational changes induced by desiccation, and reduced efficiency of recovery during rehydration [4,5]. To incorporate these effects into the model, loss parameters were introduced to account for the stability of different biomolecular species:

mRNA ($M$):prone to chemical degradation and structural cleavage; loss constant $k_1$[6].

Immature polypeptides: ($P_{\text{immature}}$):unstable due to incomplete folding and modification; susceptible to aggregation or degradation; loss constant $k_2$

Mature proteins ($P_{\text{mature}}$):relatively stable but still subject to conformational changes and partial denaturation; loss constant $k_3$[7].

Cell-wall-anchored proteins ($P_{\text{free}}$):anchoring enhances stability, but cell wall structure may be disrupted by freeze-drying stress [8]; loss constant $k_4$

Thus, the initial values for the rehydration stage were defined as:
\[
    M_0 = (1 - k_1) \cdot M
\]
\[
    P_{\text{immature},0} = (1 - k_2) \cdot P_{\text{immature}}
\]
\[
    P_{\text{mature},0} = (1 - k_3) \cdot P_{\text{mature}}
\]
\[
    P_{\text{wall},0} = (1 - k_4) \cdot P_{\text{wall}}
\]

Estimated parameter values: $k_1$ = 0.2,$k_2$ = 0.25,$k_3$ = 0.1,$k_4$ = 0.05

|Parameter          |Meaning                                                      |
|:------------------|:-------------------------------------------------------------|
|$M_0$            |Amount of fusion protein mRNA after lyophilization          |
|$P_{\text{immature}}$ |Amount of immature peptide chains after lyophilization      |
|$P_{\text{mature}}$   |Amount of mature proteins after lyophilization              |
|$P_{\text{wall}}$     |Amount of fusion proteins anchored on the cell wall after lyophilization |
|$k_1$            |Loss constant of fusion protein mRNA during lyophilization  |
|$k_2$            |Loss constant of immature peptide chains during lyophilization |
|$k_3$            |Loss constant of mature proteins during lyophilization      |
|$k_4$            |Loss constant of cell wall-anchored proteins during lyophilization |

Since no direct experimental data are available for biomolecular damage during yeast lyophilization, these constants were estimated based on molecular stability and possible damage pathways. The use of cryoprotectants further reduces molecular loss, and therefore some biomolecules may remain largely intact. Moreover, fusion proteins spontaneously released into the medium prior to freezexdrying were excluded from consideration, yielding an initial $P_{\text{free}}$ = 0 at the rehydration stage [8].

c. Rehydration Stage

During rehydration, yeast cells were assumed to rapidly regain physiological activity following reconstitution [9]. Therefore, the same system of differential equations used for the cultivation stage was applied, except that the initial values corresponded to the post-lyophilization abundances.

The gene expression during the recovery phase is shown in Fig. 3.
![Fig. 3](https://static.igem.wiki/teams/5924/drylab/resuscitation-stage.webp)

The simulation results indicate that substrate-acting strain cells restored steady-state expression levels within a short period. Throughout the rehydration phase, the number of cell-wall-anchored fusion proteins remained high, whereas the abundance of spontaneously released proteins was consistently low. These findings confirm that lyophilized yeast can be successfully reactivated for botulinum toxin detection and further validate that spontaneous release of the fusion protein does not significantly compromise detection accuracy.

^5.  Discussion:^

In this study, most parameter values were obtained from published literature. Although discrepancies with actual experimental conditions may exist due to environmental variability, the overall qualitative trends align well with biological expectations.
Two core assumptions underlie the modeling framework:

1) We assumed immediate anchoring of fusion proteins upon secretion into the periplasmic space. This simplification effectively neglects anchoring delay, which is justified given the lack of reliable anchoring rate measurements and the strong enzymatic propensity of PIR proteins to covalently attach to the cell wall.

2) We neglected the transient lag in physiological recovery during rehydration. Although in reality yeast does not immediately regain full activity, rehydrated dry yeast has been shown to recover rapidly in suitable media. Therefore, this assumption does not substantially affect the modeled dynamics.

Although several parameters were estimated rather than experimentally measured, all estimations were made within biologically reasonable ranges and do not alter the qualitative conclusions. Overall, this modeling study provides a quantitative framework for analyzing the complete workflow of substrate strain production, preservation, and application.

^6. References^

[1]:Martinić Cezar T, Lozančić M, Novačić A, Matičević A, Matijević D, Vallée B, Mrša V, Teparić R, Žunar B. Streamlining N-terminally anchored yeast surface display via structural insights into S. cerevisiae Pir proteins. Microb Cell Fact. 2023 Sep 7;22(1):174. doi: 10.1186/s12934-023-02183-2. PMID: 37679759; PMCID: PMC10483737.

[2]:Sütterlin C, Doering TL, Schimmöller F, Schröder S, Riezman H. Specific requirements for the ER to Golgi transport of GPI-anchored proteins in yeast. J Cell Sci. 1997 Nov;110 ( Pt 21):2703-14. doi: 10.1242/jcs.110.21.2703. PMID: 9427388.

[3]:Martinić Cezar T, Lozančić M, Novačić A, Matičević A, Matijević D, Vallée B, Mrša V, Teparić R, Žunar B. Streamlining N-terminally anchored yeast surface display via structural insights into S. cerevisiae Pir proteins. Microb Cell Fact. 2023 Sep 7;22(1):174. doi: 10.1186/s12934-023-02183-2. PMID: 37679759; PMCID: PMC10483737.

[4]:Iyer LK, Sacha GA, Moorthy BS, Nail SL, Topp EM. Process and Formulation Effects on Protein Structure in Lyophilized Solids Using Mass Spectrometric Methods. J Pharm Sci. 2016 May;105(5):1684-1692. doi: 10.1016/j.xphs.2016.02.033. Epub 2016 Apr 1. PMID: 27044943; PMCID: PMC4846509.

[5]:Tyagi, N.; Gidlöf, Z.; Osanlóo, D.T.; Collier, E.S.; Kadekar, S.; Ringstad, L.; Fureby, A.M.; Roos, S. The Impact of Formulation and Freeze Drying on the Properties and Performance of Freeze-Dried Limosilactobacillus reuteri R2LC. Appl. Microbiol. 2023, 3, 1370-1387. https://doi.org/10.3390/applmicrobiol3040092

[6]:Khan MDFH, Baudin F, Sudalaiyadum Perumal A, Kamen AA. Freeze-Drying of mRNA-LNPs Vaccines: A Review. Vaccines (Basel). 2025 Aug 12;13(8):853. doi: 10.3390/vaccines13080853. PMID: 40872937; PMCID: PMC12389932.

[7]:Kumar R, Kharbikar BN. Lyophilized yeast powder for adjuvant free thermostable vaccine delivery. Appl Microbiol Biotechnol. 2021 Apr;105(8):3131-3143. doi: 10.1007/s00253-021-11259-1. Epub 2021 Apr 9. PMID: 33834253; PMCID: PMC8032460.

[8]:Nowak D, Jakubczyk E. The Freeze-Drying of Foods-The Characteristic of the Process Course and the Effect of Its Parameters on the Physical Properties of Food Materials. Foods. 2020 Oct 18;9(10):1488. doi: 10.3390/foods9101488. PMID: 33080983; PMCID: PMC7603155.

[9]:Novo M, Beltran G, Rozes N, Guillamon JM, Sokol S, Leberre V, François J, Mas A. Early transcriptional response of wine yeast after rehydration: osmotic shock and metabolic activation. FEMS Yeast Res. 2007 Mar;7(2):304-16. doi: 10.1111/j.1567-1364.2006.00175.x. Epub 2006 Nov 21. PMID: 17132143.

[10]:Pelechano, V., Chávez, S., & Pérez-Ortín, J. E. (2010). A complete set of nascent transcription rates for yeast genes. PLOS ONE, 5(11), e15442. https://doi.org/10.1371/journal.pone.0015442

[11]:Y. Wang, C.L. Liu, J.D. Storey, R.J. Tibshirani, D. Herschlag, & P.O. Brown, Precision and functional specificity in mRNA decay, Proc. Natl. Acad. Sci. U.S.A. 99 (9) 5860-5865, https://doi.org/10.1073/pnas.092538799 (2002).

[12]:Pfeffer, M., Maurer, M., Köllensperger, G. et al. Modeling and measuring intracellular fluxes of secreted recombinant protein in Pichia pastoris with a novel 34S labeling procedure. Microb Cell Fact 10, 47 (2011). https://doi.org/10.1186/1475-2859-10-47

[13]:Lacroute F. RNA and protein elongation rates in Saccharomyces cerevisiae. Mol Gen Genet. 1973 Sep 27;125(4):319-27. doi: 10.1007/BF00276587. PMID: 4591362.

[14]:Sütterlin C, Doering TL, Schimmöller F, Schröder S, Riezman H. Specific requirements for the ER to Golgi transport of GPI-anchored proteins in yeast. J Cell Sci. 1997 Nov;110 ( Pt 21):2703-14. doi: 10.1242/jcs.110.21.2703. PMID: 9427388.

[15]:Yewdell JW. DRiPs solidify: progress in understanding endogenous MHC class I antigen processing. Trends Immunol. 2011 Nov;32(11):548-58. doi: 10.1016/j.it.2011.08.001. Epub 2011 Sep 29. PMID: 21962745; PMCID: PMC3200450.

[16]:Belle A, Tanay A, Bitincka L, Shamir R, O'Shea EK. Quantification of protein half-lives in the budding yeast proteome. Proc Natl Acad Sci U S A. 2006 Aug 29;103(35):13004-9. doi: 10.1073/pnas.0605420103. Epub 2006 Aug 17. PMID: 16916930; PMCID: PMC1550773.

[17]:Hirschberg K, Miller CM, Ellenberg J, Presley JF, Siggia ED, Phair RD, Lippincott-Schwartz J. Kinetic analysis of secretory protein traffic and characterization of golgi to plasma membrane transport intermediates in living cells. J Cell Biol. 1998 Dec 14;143(6):1485-503. doi: 10.1083/jcb.143.6.1485. PMID: 9852146; PMCID: PMC2132993.

[18]:Mahbuba H. Meem, Paul J. Cullen, The impact of protein glycosylation on Flo11-dependent adherence in Saccharomyces cerevisiae, FEMS Yeast Research, Volume 12, Issue 7, November 2012, Pages 809–818, https://doi.org/10.1111/j.1567-1364.2012.00832.x

# The entire process simulation of the project: the corresponding relationship between botulinum toxin concentration and signal strength

^1. Motivation^

In this project, we aim to design a botulinum neurotoxin A (BoNT/A) biosensor based on the chassis organism *Saccharomyces cerevisiae*. Utilizing BoNT/A's cleavage activity on SNAP25, the cascade amplification of the MATα mating type's MAPK signaling pathway in response to α-factor, and the downstream product MjDOD catalyzing a reaction to produce the final colored product betaxanthin (yellow), we hope to achieve a sensitive and intuitive portable detection method for BoNT/A.

However, the following key issues cannot be resolved by experiments alone:

a.Since our laboratory cannot directly use BoNT/A as an experimental substance, how will we ultimately determine the relationship between BoNT/A concentration and the final signal?

b.Due to the inherent complexity of biological systems and the relative scarcity of high-quality data, there is significant uncertainty in model parameters. Although parameter optimization cannot be performed before obtaining rigorous experimental data, how can we understand the impact of different parameters on the prediction results to minimize errors?

c.The Limit of Detection (LOD) and Limit of Quantification (LOQ) are crucial attributes of a biosensor and key differentiators for its precision. However, since our team cannot obtain BoNT/A, how can we determine the designed LOD and LOQ?

To address these problems, as a substitute for real experiments, to derive the complete BoNT/A–signal relationship, calculate the LOD and LOQ, attempt to compare them with the detection thresholds of the human eye and micro-spectrophotometers, and better guide subsequent design to avoid blind screening of experimental conditions, we have established a mathematical model.

^2.Overview^

![](https://static.igem.wiki/teams/5924/drylab/graph-demo.webp)

Our model is primarily divided into three parts: Plotting the BoNT/A–signal dose-response curve, Stochastic process simulation and calculation of LOD/LOQ, and Parameter sensitivity analysis. The latter two are analyses of the former's predictions and are described in detail elsewhere. The former, "BoNT/A–signal dose-response curve plotting," directly involves the entire pathway and is relatively complex; therefore, we provide a detailed analysis here:

^a.^ Two types of *S. cerevisiae*, the upstream strain (Substrate-acting strain) and the downstream strain (Signal-transferring Strain), are co-cultured in liquid medium, competing until they reach equilibrium and the total concentration stabilizes. The upstream strain expresses a fusion protein of SNAP25 C-terminally fused to α-factor on its surface, while the downstream strain has the natural α-factor receptor Ste2 on its surface.

^b.^ We engineered the MAPK signaling pathway in the upstream strain: knocking out the *ste12*, *gal4*, *gal80*, and *far1* genes (*ste12Δ, gal4Δ, gal80Δ, far1Δ*), and expressing a synthetic fusion protein, sTF, inside the cell.

^c.^ Binding of α-factor to the Ste2 receptor activates this signaling pathway. The downstream Fus3 protein in the pathway can interact with the exogenous fusion protein sTF. sTF is a fusion protein with a *Ste12PRD* at one end and a *Gal4DBD* at the other end. Fus3 protein activates sTF, which then, via *Gal4DBD*, activates the *pGAL1* promoter in the cell. We knocked out the endogenous *GAL1* gene and integrated the gene for *DODA* (4,5-dopa dioxygenase) from *Mirabilis jalapa* (*MjDOD*) downstream of the *pGAL1* promoter into the yeast nuclear genome. We expect that after the Fus3 protein activates the *pGAL1* promoter, it will activate Gal4DBD and successfully express the MjDOD protein.

^d.^ At the start of detection, we add the sample to be tested, which may contain trace amounts of BoNT/A, to the culture medium in the detection kit. In a theoretical laboratory validation, we would add a quantified amount of BoNT/A to the medium. It exerts its protease activity, cleaving the SNAP25 part of the fusion protein and releasing the α-factor. The α-factor diffuses and binds to the Ste2 receptors on the surface of MATα yeast, transmitting the signal downstream via the MAPK signaling pathway, initiating downstream signaling and metabolic pathways, leading to the synthesis of MjDOD. Simultaneously, the α-factor and Ste2 receptors are internalized and degraded.

^e.^ Concurrently, we add a certain amount of L-DOPA to the medium. MjDOD, acting as a dopa dioxygenase, catalyzes the conversion of intracellular L-DOPA into betalamic acid. Betalamic acid then spontaneously reacts with various amino acids in the medium to generate different betaxanthins.

^f.^ One method of laboratory validation involves using a spectrophotometer to measure absorbance and calculating the betaxanthin concentration according to the Beer-Lambert law. By detecting the dose-response relationship between the added BoNT/A concentration and the final signal at different time points, the performance of the biosensor is ultimately validated.

^3. Modeling Objective & Approach^

We constructed a mathematical model primarily based on Delay Differential Equations (DDEs) to simulate the entire process from BoNT/A cleavage of the fusion protein to the final light absorption by betaxanthins.

We introduced delay times for steps that clearly take significant time, such as signal transduction and mRNA nuclear retention. We considered background noise and measurement errors, and introduced Gaussian random noise at key steps for Monte Carlo simulations, striving to enhance the simulation's credibility under the challenging condition of limited reference data.

Our model aims to achieve at least the following objectives:

a. Derive the relationship curve between the betaxanthin signal intensity and experimental time, initial BoNT/A concentration, and the initial ratio of the two yeast strains.

b. Calculate the LOD and LOQ of this biosensor and provide the variation of LOD and LOQ under different conditions.

After establishing the basic model framework, we initiated a new "Design-Build-Test-Learn" cycle to improve the model's predictive accuracy, achieving the following objectives:

a. Evaluate the optimal yeast ratio at the beginning of the plateau phase to predict the best initial yeast ratio and optimize culture conditions.

b. Evaluate the parameter sensitivity of the prediction results to provide a basis for further parameter adjustment.

c.To build  an easily-used interactive computational platform for modeling, simulating, and analyzing the limit of detection in biological signaling pathways without writing any code.

We used ^optimization algorithms^ to scan different concentration ratios of the upstream and downstream strains, striving to find the initial ratio that minimizes the LOD and maximizes the signal. We designed our own algorithm to evaluate the sensitivity of the model predictions to each parameter individually and used the ^Sobol index method^ to perform a global sensitivity analysis of the model parameters.

After establishing the model, we concurrently conducted further research. Noting the similarities between our algorithm and existing modeling software like ^COPASI^ , but considering COPASI's shortcomings such as high learning curve, relative complexity, poor cross-platform compatibility, limited customizability, and complex visualization and interaction experience, and also to promote the dissemination of results within the team and to external audiences, and to facilitate quick and intuitive presentation of product effects to users, we also endeavored to build ^a web platform^. This platform not only serves as a display for the model but, due to its ability to flexibly add/remove nodes and edges (interactions), holds promise to fill the gap as a learning and rapid experimental platform for future teams modeling biological systems.

^4. Assumptions^

^Homogeneity^: 

Except where replaced by delay times, the intracellular and extracellular environments are well-mixed systems without spatial gradients.

^Stability^: 

Within the simulation time scale, the concentration and enzymatic activity of BoNT/A remain constant. Yeast cell growth has entered the stationary phase, the total concentration remains constant, and various states are relatively stable.

^Partial Quasi-Steady State^: 

The binding of α-factor and the Ste2 receptor is assumed to be in a quasi-steady state, maintaining binding and dissociation equilibrium instantaneously.

^Limited Determinism/Uncertainty^: 

Biological systems are highly complex, and the parameters used for modeling themselves have significant uncertainty, thus it is necessary to introduce randomness. However, due to the excessive uncertainty inherent in biological systems, performing uncertainty modeling for every process is too workload-intensive. To minimize model complexity, we only perform uncertainty modeling for the most critical processes like signal transduction, transcription, and translation, and make deterministic assumptions for other processes with lower uncertainty.

^Limited Spatial-Temporal Delay^:

Since biological reactions always occur asynchronously and potentially in different locations, but to simplify the model as much as possible, we assume that significant time delay effects due to diffusion etc., exist only for key processes like α-factor diffusion, signal transduction, and mRNA nuclear retention/export. For other processes, reactions are assumed to be synchronous. Delay effects in other processes are neglected.

^Stable System Assumption^: 

Effects of environmental fluctuations such as temperature and pH, other than changes in the concentrations of some medium components, on the model are not considered.

^Neglect of Secondary Processes and Factors^: 

Enzymatic reactions follow ^Michaelis-Menten kinetics^ , without cooperative effects or product inhibition.

The absorbance of betalamic acid, L-DOPA, etc., at λ = 480nm is neglected. (This assumption is acceptable based on our experimental data; the absorbance of blank medium with experimental concentrations of L-DOPA is very low (Abs ≈ 0.01), and betalamic acid has low absorbance at λ = 480nm. In the absence of more detailed real data, this assumption is acceptable).

Due to lack of data, differences in the binding of various amino acids to betalamic acid are ignored, and amino acids are treated as an "average amino acid." 

Due to the strong affinity between L-DOPA and MjDOD, which greatly exceeds non-catalytic collisions between L-DOPA and betalamic acid, the binding of L-DOPA and betalamic acid is neglected. The binding of substances other than amino acids with betalamic acid is neglected.

^5. Equations^

| Variable                               | Description                                                                   | Unit     |
| :------------------------------------- | :---------------------------------------------------------------------------- | :------- |
| $c_{\alpha{-factor-free}}$      | Free α-factor concentration                                                   | $M$    |
| $c_{{snap25}}$                  | SNAP25 concentration                                                          | $M$    |
| $n_{{ste2-receptor}}$           | Number of Ste2 receptors per downstream cell                                  | $mol/cell$ |
| $c_{{MjDOD-mRNA}}$              | Concentration of MjDOD mRNA in the system                                     | $M$    |
| $c_{{MjDOD-protein}}$           | Concentration of MjDOD protein inside downstream cells in the system          | $M$    |
| $c_{{L-DOPA}}$                  | L-DOPA concentration                                                          | $M$    |
| $c_{{aminoacids}}$              | Amino acid concentration                                                      | $M$    |
| $u_{{pGAL1}}$                   | Activation ratio of pGAL1 promoter relative to maximum activation intensity   | \-       |
| $c_{{betalamicacid}}$           | Betalamic acid concentration                                                  | $M$    |
| $c_{{betaxanthins}}$            | Betaxanthin concentration                                                     | $M$    |
| $A$                                  | System absorbance intensity                                                   | \-       |

^a. $BoNT/A\ Cleavage \to GAL1\ Promoter\ Activation$^

\[
\begin{cases}
\begin{cases}
inflow\_rate\_alphafactor  = \frac{k_{cat_{bonta}}*c_{snap25}}{K_{m_{bonta}}+c_{snap25}}*c_{bonta}*c_{snap25}\\
outflow\_rate\_alphafactor = - c_{\alpha-factor-free}*\frac{\ln2}{t_{0.5-\alpha factor}}\\
\frac{dc_{\alpha-factor-free}}{dt} = inflow\_rate\_alphafactor \ +\  outflow\_rate\_alphafactor\\
\frac{dc_{snap25}}{dt} = -\frac{k_{cat_{bonta}}*c_{snap25}}{K_{m_{bonta}}+c_{snap25}}*c_{bonta}*c_{snap25}\\
\frac{c_{ste2-receptor}}{dt} = - c_{\alpha-factor-free}*\frac{\ln2}{t_{0.5-\alpha factor}}\\
\theta = \frac{c_{\alpha-factor-free}}{EC50_{ste2}+c_{\alpha-factor-free}}\\
u_{pGAL1}(t) = \theta(t - t_{delay-transduction})\\
\end{cases}\\,
\begin{cases}
c_{snap25}(0) = c_{snap25-initial}\\
c_{ste2-receptor}(0) = c_{ste2-receptor-initial} \\
c_{\alpha-factor-free}(0) = c_{\alpha-factor-free-initial}\\

\end{cases}
\]

\[
\begin{cases}
c_{snap25-initial} = \frac{\frac{ration\_cell}{ration\_cell+1}*c_{cell}*n_{snap25}}{N_{A}}\\
c_{ste2-receptor-initial}=\frac{\frac{1}{ration\_cell+1}*c_{cell}*n_{ste2-receptor}}{N_{A}}\\
\end{cases}
\]

BoNT/A cleaves the SNAP25 part of the fusion protein on the surface of MATA-type yeast, releasing α-factor. The α-factor diffuses to the vicinity of MATα yeast, binds to the Ste2 receptor, and simultaneously, according to the concentration of α-factor bound to Ste2 receptors, influences the internalization and degradation of α-factor and receptors.

^b. $GAL1\ Promoter\ Driven\ mRNA\ Transcription \to MjDOD\ Protein\ Synthesis\ and\ Degradation$^

\[
\begin{cases}
k_{transcription-pGAL1} = k_{transcription-pGAL1-max}*u_{pGAL1}\\
inflow\_rate\_mjdod\_mrna= k_{transcription-pGAL1}\\
\frac{dc_{mjdod-mrna}}{dt}  = inflow\_rate\_mjdod\_mrna(t_{delay-transduction}<t<t_{delay-transduction}+t_{delay-transcription-to-translation})
\end{cases}
\]

\[
c_{mjdod-mrna}(0) = c_{mjdod-mrna-initial}
\]

The transcription rate of the *pGAL1* promoter is calculated proportionally to its maximum transcription rate. Downstream MjDOD mRNA is transcribed at this rate.

\[
\begin{cases}
inflow\_rate\_mjdod\_protein(t) = k\_translation\_mjdod*c_{mjdod-mrna}\\
outflow\_rate\_mjdod\_protein(t) = - c_{mjdod-protein}*\frac{\ln2}{t_{0.5-mjdod-protein}}\\
\frac{dc_{mjdod-protein}}{dt} = inflow\_rate\_mjdod\_protein\ + \ outflow\_rate\_mjdod\_protein \\
outflow\_rate\_mjdod\_mrna = - c_{mjdod-mrna}*\frac{\ln2}{t_{0.5-mjdod-mrna}}\\
\frac{dc_{mjdod-mrna}}{dt} = inflow\_rate\_mjdod\_mrna \ +\ outflow\_rate\_mjdod\_mrna
\end{cases}
\]

After MjDOD mRNA enters the cytoplasm, it is translated to produce MjDOD protein. Both mRNA and protein begin degradation simultaneously.

^c. $MjDOD\ Catalysis \to Betaxanthin\ Color\ Development$^

\[
\begin{cases}\frac{dc_{ldopa}}{dt} = - \frac{k_{cat_{mjdod}}}{K_{m_{mjdod}}+c_{ldopa}}*c_{ldopa}*c_{mjdod-protein}
\\inflow\_rate\_betalamic\_acid = \frac{k_{cat_{mjdod}}}{K_{m_{mjdod}}+c_{ldopa}}*c_{ldopa}*c_{mjdod-protein}
\\degradation\_rate\_betaxanthins= - c_{betaxanthins}*\frac{\ln2}{t_{0.5-betaxanthins}}
\\\frac{dc_{aminoacids}}{dt} = outflow\_rate\_betalamic\_acid =  - k_{betalamicacid\_to\_betaxanthins}*c_{aminoacids}*c_{betalamic-acid}
\\\frac{dc_{betalamic-acid}}{dt} = inflow\_rate\_betalamic\_acid\ +\ outflow\_rate\_betalamic\_acid
\\\frac{dc_{betaxanthins}}{dt} = k_{betalamicacid\_to\_betaxanthins}*c_{aminoacids}*c_{betalamic-acid}\ +\ degradation\_rate\_betaxanthins
\end{cases}
\]

MjDOD catalyzes the conversion of L-DOPA to betalamic acid; betalamic acid then combines with amino acids inside the cell to generate (mixed) betaxanthins. Betaxanthins produce color and simultaneously undergo degradation.

^d.$Monte\ Carlo\ Method\ for\ Noise\ Simulation$^

\[
\begin{cases}
for\ final\  signal:
\\
\sigma_{final-signal} = A_{final}(t,c_{\alpha-factor}) \cdot CV_{final-signal}
\\
for \ blank \ signal:
\\
\sigma_{blank-signal} = A_{blank}(t,c_{\alpha-factor}) \cdot CV_{blank-signal}
\\assuming :u_1, u_2 are random number in (0,1) ,then:\\z_0 = \sqrt{-2 \ln u_1} \cdot \cos(2\pi u_2)
\\{noise} = z_0 \cdot \sigma + \mu
\end{cases}
\]

The Box-Muller method is used to simulate the random noise distribution of the final signal and background noise. Assuming it follows a normal distribution based on the characteristics of the biological system.

^e. $Limit \ of\  Detection\  (LOD) \ \& \  Limit\  of\  Quantification\  (LOQ) \ Calculation$^

\[
\begin{cases}
Threshold_{LOD}=μ_{blank} +3⋅σ_{blank}
\SNR _{LOD}= \frac{S(c_{LOD})−μ_{blank}}{σ_{blank}}=3)
\\Threshold_{LOQ}=μ_{blank}+10⋅σ_{blank}
\SNR _{LOQ}= \frac{S(c_{LOQ})−μ_{blank}}{σ_{blank}}=10)
\\c_{{LOD}} = f^{-1}\left( {Threshold}_{{LOD}} \right)
\\c_{{LOQ}} = f^{-1}\left( {Threshold}_{{LOQ}} \right)
\end{cases}
\]

Using local interpolation, obtain the inverse mapping $f^{-1}$ from BoNT/A concentration to signal, then calculate LOD and LOQ according to the IUPAC recommended method above.

^f. $Single\ Parameter\ Sensitivity\ Analysis$^

\[
\begin{cases}
S= \frac{ \  Relative\  Variation_{final-signal}}{\  Relative\ Change_{signal-parameter}}= \frac{\frac{Δy_{avg}}{y_{0}}}{\frac{Δp_{avg}}{p_{0}}}\\
Δy_{avg} = \frac{y_{+}-y{-}}{2}\\
Δp_{avg} = \frac{p_{+}-p{-}}{2}\\
\end{cases}
\]

The Sensitivity Coefficient $S$ represents the relative strength of the relationship between parameter change and output change. $S > 0$ indicates a positive correlation; $S < 0$ indicates a negative correlation; the larger $|S|$ is, the more significant the parameter's influence on the output.

^g. $Global\  Sensitivity\  Analysis$^

The Sobol method, based on variance decomposition, quantifies the contribution of parameters and their interactions to the output variation. The formulas are as follows:

\[
\begin{cases}
V(Y) = \sum_{i=1}^n V_i + \sum_{1 \leq i<j \leq n} V_{ij} + \sum_{1 \leq i<j<k \leq n} V_{ijk} + \dots + V_{1,2,\dots,n} \\
S_1^{(i)} = \frac{V_i}{V(Y)}, \quad ST^{(i)} = \frac{V_i + \sum_{j \neq i} V_{ij} + \dots + V_{1,2,\dots,n}}{V(Y)}
\end{cases}
\]

$V(Y)$ represents the total variance of the output variable; $V_i$ represents the variance caused by parameter $i$ alone; $V_{ij}$ represents the variance caused by the interaction between parameters $i$ and $j$. The first-order index $S_{1_i}$ indicates the proportion of output variation contributed by parameter $i$ alone. Similarly, the total effect index $ST_i$ represents the total proportion of output variation contributed by parameter $i$ and its interactions with all other parameters.

Both $S_1$ and $ST$ are non-negative; a larger value indicates a more significant influence of the parameter (or including its interactions) on the output. $S_{T_i} - S_{1_i}$ reflects the strength of the interaction between parameter $i$ and other parameters.

^f. Initial Values of Variables^

| Initial Value of Variables                                | Description                                                                                             | Value               | Reference                       |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | ------------------- | ------------------------------- |
| $c_{α-factor-initial}$           | Initial α-factor concentration                                                                          | 0                   | Project design                  |
| $n_{MjDOD-gene}$                         | Copy number of MjDOD gene per downstream cell                                                           | $1{mol/cell}$ | Project design                  |
| $c_{MjDOD-mRNA-initial} $         | Initial MjDOD mRNA concentration *More precise parameters are available here; updating now is recommended to reduce noise and improve detection limit* | $0.1{mol/cell}$ | (V. Iyer et. al.,1996)[1]|
| $c_{LDOPA-initial}  $ | Initial L-DOPA concentration | $10^{-3} M$ | Project design |
| $c_{aminoacids-initial} $| Initial intracellular amino acid concentration  | $9.88 \times 10^{-3} M$ | Culture medium instructions |
| $A_{initial} $ | Initial absorbance intensity | 0 | Estimated |
|  $c_{betalamic-acid-initial}$  | Initial betalamic acid concentration  | 0  |  Estimated |
|  $ c_{betaxathins-initial}$    | Initial betaxanthin concentration | 0 |  Estimated  |
|  $c_{MjDOD-initial}$ | Initial MjDOD protein concentration|0 |  Estimated  |
| $n_{ste2-receptor-initial} $   | Initial number of Ste2 receptors per upstream cell surface | $8000\ {mol/cell}$ | (Jenness, D. D. et. al.,1986)[2] |
|  $c_{cell}$ | Total cell concentration of both strains | $2 \times 10^{11}\ {cell/L} $| [3]     |
|ration_cell                             | Ratio of upstream to downstream cell numbers | \ | Manually set  |
|  $c_{bonta}$  | BoNT/A concentration | \                   | Measurement target              |
| $u_{pGAL1-initial}$ | Ratio of uninduced pGAL1 promoter transcription rate to fully induced pGAL1 transcription rate          | 0.1\%           | [21]|
|  $n_{snap25-initial} $  | Initial number of SNAP25 per upstream cell surface | $3256 ± 887 mol/cell$ | Saccharomyces genome database[4] |

^6.Parameters^

| Parameter                                | Description                                                  | (Adopted) Value                                              | Reference                                                    |
| ---------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| $EC50_{Ste2}$                            | EC50 for downstream mRNA transcription rate upon α-factor binding to Ste2 receptor | $8.3 \times 10^{-11} M $                                     | (Fan, C. et. al.,2024)[5]                                    |
| $K_{m_{mjdod}} $                         | K_m for MjDOD catalyzing L-DOPA to betalamic acid            | $1.6848 \times 10^{-4} M$                                    | (Chou YC et. al.,2019)[6]                                    |
| $ k_{cat_{mjdod}} $                      | k_{cat} for MjDOD catalyzing L-DOPA to betalamic acid        | $2.12*10^{-2}s^{-1}$                                         | (Chou YC et. al.,2019)[6]                                    |
| $K_{m_{bonta}} $                         | K_m for BoNT/A protease activity (1h)                        | $ 1.89±0.146*10^{-4}M$                                       | (Mizanur RM et. al.,2014)[7]                                 |
| $k_{cat_{bonta}}$                        | k_{cat} for BoNT/A protease activity (1h)                    | $134 ±8.75s^{-1}$                                            | (Mizanur RM et. al.,2014)[7]                                 |
| $M_{bonta}$                              | Molecular weight of BoNT/A core (heavy + light chains)       | $149425.84Da(149\ kDa)$                                      | Expasy protparam calculation[8]                              |
| NA $N_A$                                 | Avogadro's number                                            | $6.022*10^{23}mol^{-1}$                                      | \                                                            |
| $OD_{600}$                               | Cell concentration corresponding to 1 OD600 unit at 600nm    | $1.4*10^{7}cell/mL (1.4*10^{10}cell/L)$                      | [9]                                                          |
| $k\_transcription\_pGAL1\_blank	$     | Leakage transcription rate of GAL1 promoter                  | $0.02 per/min$                                               | (Hsu, C. et al.,2012)[10]                                    |
| $t_{0.5-\alpha factor}$                  | Half-life of bound α-factor degradation (Ste2 receptor internalization) | $4.5min (270s)$                                              | (T. Yi et. al.,2013)[11]                                     |
| $\epsilon_{betaxanthins}$                | Maximum absorption peak and average molar extinction coefficient of various betaxanthins | $4.8*10^{6}M^{-1}m^{-1}$                                     | (Girod, PA. et. al.,1991)[12]                                |
| $c_{ss-GAL1mRNA-gal(+)}$                 | GAL1 mRNA number per yeast cell with galactose addition      | $33±5 per/cell$                                              | (V. Iyer et. al.,1996)[1]                                    |
| $halflife_{GAL1mRNA}$                    | GAL1 mRNA degradation half-life                              | $8.1min$                                                     | (Hsu, C. et al.,2012)[10]                                    |
| $k_{transciption-pGAL1-max}$             | Maximum induced transcription rate of pGAL1 promoter         | 2.82\ {molecules/(cell·min)}                                 | estimated(V. Iyer et. al.,1996, Hsu, C. et al.,2012)^{1} [1][10] |
| $k\_transcription\_pGAL1\_blank	$     | Leakage transcription rate of pGAL1 promoter                 | 0.02\ {molecules/(cell·min)}                                 | (Hsu, C. et al.,2012)[10]                                    |
| $rate_{ribo-elongation} $                | Average mRNA translation elongation rate in yeast            | $9 codon*script^{-1}s^{-1}$                                  | (Gao, X. et. al.,2025)[13]                                   |
| $length_{mjdod-protein}$                 | Amino acid count of MjDOD protein                            | 267                                                          | Uniprot database[14]                                         |
| $density_{ribo-per-kb}$                  | Average ribosome density per kb of yeast mRNA                | $0.64\ ribo/100b$                                            | (Steinbrecht, D. et. al.,2024)[15]                           |
| $k\_translation\_mjdod$                  | MjDOD translation rate                                       | $10.368\  protein/(mrna*min)$                                | estimated(Gao, X. et. al.,2025, Uniprot)^{2} [13][14]        |
| $k_{betalamicacid\_to\_betaxanthins}$    | Rate constant for betalamic acid reaction with amino acids to form betaxanthins | $189h^{-1}M^{-1}\star\\(3.15min^{-1}M^{-1},0.0525s^{-1}M^{-1})$ | estimated(Contreras-Llano LE et. al.,2019)^{3} [16]          |
| $optical\_distance$                      | Optical path length for spectrophotometric measurement       | $0.01m$                                                      | Experimental setting                                         |
| $t_{0.5-mjdod-mrna}$                     | MjDOD mRNA degradation half-life                             | $4.8min(288s)$                                               | (Chan, L. Y. et. al.,2018)[17]                               |
| $t_{0.5-mjdod}$                          | MjDOD protein degradation half-life                          | $20h (1200min ,72000s)$                                      | estimated(Expasy protparam)                                  |
| $t_{delay-transduction}$                 | Time delay for signal transduction process                   | $22.5min(1350s)$                                             | estimated(Aymoz, D. et. al.,2018)^{4} [18]                   |
| $t_{delay-transcription-to-translation}$ | Time delay from mRNA transcription completion to nuclear export and ribosome binding | 78\ {min} (4680\ {s}) *                                      | estimated(Steinbrecht, D. et. al.,2024)^{5} [19]             |
| $t_{0.50-betaxanthins}$                  | Betaxanthin degradation half-life                            | $297.6h\\(1.79*10^{4}min , 1.07*10^{6}s)$                    | estimated(Cruz, S et. al.,2014)^{6} [20]                     |
| $CV_{final-signal}$                      | Coefficient of variation for final signal intensity due to all factors | $40\%$                                                       | estimated                                                    |
| $CV_{final-signal}$                      | Coefficient of variation for background signal intensity due to all factors | $40\%$                                                       | estimated                                                    |
| $CV_{transduction} $       | Coefficient of variation for promoter activation intensity due to MAPK signal transduction | $30\%$                                                       | estimated                                                    |
| $CV_{transcription}$ | Coefficient of variation for mRNA transcription rate         | $30\%$                                                       | estimated                                                    |
| $CV_{translation}$    | Coefficient of variation for protein translation rate        | $20\%$                                                       | estimated                                                    |
| $CV_{measurement}$                       | Coefficient of variation for spectrophotometric measurement  | 0.002 A (1.0 mm path length, ~0.33% at 0.60A) or 1% CV, whichever is larger. Here 1% is chosen as a conservative estimate | [11]                                                         |

(1):

Calculated as:

\[
k_{transcription-pGAL1-max} = \frac{\ln 2 \times c_{ss-GAL1mRNA}}{halflife_{GAL1mRNA}} = \\ \frac{0.693 \times 33\ {molecules/cell}}{8.1\ {min}} = 2.82\ {molecules/(cell*min)}
\]

(2):

Since the total reaction time required is long (≈8h), it can be assumed that mRNA transcription initiation and elongation reach a quasi-steady state after a short period. Furthermore, the protein production rate per single mRNA is primarily determined by the mRNA-specific initiation rate, which dictates the frequency of new ribosomes starting translation. Once steady state is reached, the rate of ribosome completion equals the initiation rate.

Calculation of initiation rate:

The rate of protein production during translation is primarily determined by the process of ribosomes moving along the mRNA. Therefore:

\[
n_{ribo-per-script}=\frac{length_{mjdod-mrna-ORF}}{100}*{density_{ribo-per-kb}}\\
t_{initial-transaltion-gap-per protein}= \frac{length_{mjdod-mrna-ORF}}{3*rate_{ribo-elongation}*n_{ribo-per-script}}\\
rate_{mjdod-protein-synthesis-per-mrna} =\frac{60}{t_{initial-transaltion-gap-per protein}}
\]

$n_{ribo-per-script}$ refers to the number of ribosomes bound per open reading frame (ORF), determined by ORF length and a relatively constant ribosome "density". Let the ORF length be $length_{mjdod-mrna-ORF}$. Ribosome density is $density_{ribo-per-kb}$.

$rate_{ribo-elongation} = 9\ {codons/(transcript·s)}$[13]

ORF length:

$length_{mjdod-mrna-ORF} = 3 \cdot length_{mjdod-protein} + 3 = 267 \cdot 3 + 3 = 804\ {nt}$(Uniprot Dataset)

Average ribosome density:

$density_{ribo-per-kb} = 0.64\ {ribosomes/100 nt}$

^Annotations: Ribosome density is significantly negatively correlated with ORF length. Short ORFs (<400 nt): average density 1.2 / 100 nt; Long ORFs (>3600 nt): average density 0.14 / 100 nt.^[15]

Based on the above, the translation rate is calculated as:

\[
k\_translation\_mjdod = \frac{60}{\frac{804}{3*9*\frac{804*0.64}{100}}} =10.368\ protein/(mrna*s)
\]

(3):
    
In the severe lack of relevant data, we can only assume that the binding of L-DOPA and betalamic acid is similar to that of general amino acids with betalamic acid. The following facts and conclusions support our assumption:
    
L-DOPA, like general amino acids, is an amino-containing molecule. The reaction mechanism with betalamic acid involves spontaneous condensation to form an imine structure (Schiff base), ultimately generating betaxanthins.

(4):

![](https://static.igem.wiki/teams/5924/drylab/max-expression-response-time-fig.webp)

early: $15\ {min}\\ late: 30\ {min}$

Based on the relationship between average expression output and response time in the graph, the average response time, i.e., the signal transduction delay time, can be estimated as:$ t_{delay-transduction} ≈ 22.5\ {min} = 1350\ {s}$

(5):

"Nuclear retention half-life" means: the time required for 50% of the nuclear mature mRNA population to complete their nuclear residence. The decrease in the number of mature RNAs in the nucleus is mainly due to export and degradation. For mRNAs of most cytoplasmic proteins, nuclear degradation is almost negligible. Therefore, the "nuclear retention half-life" can be considered as the median time for mRNA export from the nucleus. Thus, it is estimated:

$t_{half\_out} = t_{delay-transcription-to-translation} = 78\ {min} = 4680\ {s}$

(6):

In the CIELab color space, the $b*$ value specifically represents the yellow-blue color component. A positive $b*$ indicates a color similar to yellow, with a larger value indicating a more intense yellow hue; a negative $b*$ indicates a color to similar to blue, with a smaller value (larger absolute value) indicating a more intense blue hue.

Coincidentally, the absorption wavelength of various betaxanthins is about 480nm, corresponding to blue-green light, while the reflected light is just yellow. Furthermore, the literature reports processing each betaxanthin separately, implying single substances. Experimental data also shows that $L*$ and $a*$ change very little before and after treatment, while only $b*$ changes significantly. In summary, we can accept using the change in $b*$ chroma to approximate represent the concentration of betaxanthins.

Due to the special color change of L-Phe-Bx, we carefully examined the L*, a*, and b* of the other three betaxanthins before storage (BS) and after storage (AS) for 8 days (L-Leu-Bx for 4 days at $ pH$=4) at $pH$=4.0 and $pH$=6.0, respectively.



![betaxanthins-lab](https://static.igem.wiki/teams/5924/drylab/betaxanthins-lab-as-bstable.webp)

| Betaxanthins | pH  | Storage time/day | b^{*}(BS) | Δb^{*} | \frac{Δb*}{b*} \times 100\% | t_{0.5}/day |
| ------------ | --- | ---------------- | ------------- | ---------- | ------------------------------- | --------------- |
| L-Ala-Bx     | 4   | 8                | 37.2          | -20.5      | -55.1                           | 6.9             |
|              | 6   | 8                | 50.9          | -7.9       | -15.5                           | 33.2            |
| L-Leu-Bx     | 4   | 4                | 35.4          | -25.4      | -71.8                           | 2.2             |
|              | 6   | 8                | 42.2          | -15.3      | -36.3                           | 12.4            |
| L-Pro-Bx     | 4   | 8                | 43.9          | -18.6      | -42.4                           | 10.1            |
|              | 6   | 8                | 47.4          | 0.4        | 0.01                            | 550             |

Observation readily shows that the half-lives of various betaxanthins are long, far exceeding the test time (8h). Therefore, the degradation of betaxanthins during the experiment can be ignored. Since our product is designed to operate in a near-neutral pH environment (pH≈6), it is reasonable to use the degradation of L-Leu-Bx under this condition to represent the overall very low degradation rate:

\[
t_{0.5 - betaxanthins} = 12.4\ {day} = 297.6\ h = 1.79 \times 10^{4}\ {min} = 1.07 \times 10^{6}\ {s}
\]

To verify the accuracy of the above estimate, we performed the following validation:

  

![betaxanthinsdegradation-hight-condition](https://static.igem.wiki/teams/5924/drylab/betaxanthinsdegradation-hight-condition.webp)

![](https://static.igem.wiki/teams/5924/drylab/betaxanthins-degradation-condition.webp)

It is easy to observe that the degradation half-life of betaxanthins at $pH$≈6.48, $T$=50.0°C is 17.41h. Based on the relationship between half-life and temperature, the half-life at room temperature (about 25°C) at the same pH will be much longer than 17h, supporting that within the modeling error tolerance, we can indeed ignore their degradation (when $t<8h$) or set it to a relatively long time.

^7. Algorithm, Implementation  & Reproducibility^

To implement the model, we comprehensively used the following algorithms to make the model more accurate and robust:

+ ^a. Innovatively introduced graph theory and graph algorithms for representation and computation.^

    We modeled based on a graph model. Nodes in the graph correspond to substances, and edges correspond to interactions between two adjacent nodes. We used Dijkstra's algorithm to calculate the sum of appropriate specific path delay times for each node to determine their earliest signal output time.

+ ^b.Adopted a "flux method" for modeling node concentration changes.^

    For each node, its own concentration (intensity) is stored. For general nodes, the input and output rates of the node's concentration (intensity) can be determined based on the Delay Differential Equation (DDE) equations of upstream and downstream edges. Likening a node to a reservoir, the node's input and output rates can be compared to the inflow rate of the inlet pipe and the outflow rate of the outlet pipe, respectively.

+ ^c.Edge process modeling primarily based on Delay Differential Equations (DDE).^

    Except for relationships between two nodes that can be expressed using elementary function expressions, for general edges connecting two nodes, they are defined by Delay Differential Equations (DDE). However, since the delay time for most processes is zero, they can also appear as ODE equations. The specific implication of the delay effect is: for an interaction process A \to B, the input/output rate calculated from the concentrations (intensities) of A and B at time t takes effect for A at time t, but for B, it only takes effect at time $t$ +$ t_{delay-time}$.

+ ^d. System dynamic simulation based on global iteration using Euler's method.^

    Using Euler's method, with a fixed time interval $Δt$, the DDE equations are converted into difference equations for iteration: $f(t+Δt) = f(t) + (rate_{inflow} - rate_{outflow}) * Δt$. We separately handled processes represented by elementary functions in the upstream (or downstream) interactions of nodes (setting their inflow or outflow rates to 0) in computation but handled them uniformly in iteration. This allows predicting the situation at every Δt interval in the future.

+ ^e.Stochastic process simulation based on the Monte Carlo method.^

    For each stochastic process, normally distributed random data were generated 100 times using the Box-Muller method, and the signal system behavior was simulated. Acceptable approximate results were obtained by approximating with probability statistics.

+ ^f. Limit of Detection (LOD) and Limit of Quantification (LOQ) calculation using local interpolation.^

    Since 3\sigma_{blank}, 10\sigma_{blank}, and \mu_{blank} do not differ greatly numerically, linear interpolation was performed near the blank signal intensity to solve for the BoNT/A concentrations corresponding approximately to \mu_{blank} + 3\sigma_{blank} and \mu_{blank} + 10\sigma_{blank}. Within the allowable error range, these are sufficient as the LOD and LOQ.

+ ^g.Single-variable sensitivity analysis.^

    For each parameter to be analyzed, we set a perturbation of ±\delta = 10\%, obtaining perturbed parameter sets. Using each perturbed parameter with other unperturbed parameters, the DDE equations were solved to obtain the model output after perturbing a single parameter, which was compared with the unperturbed output.

    The sensitivity coefficient S was calculated using the unperturbed and perturbed output results. According to the formula for S, S represents the relative strength of the absolute value's influence. A larger absolute value of S means the change in S has the greatest impact on the final signal; the sign of S indicates the direction of the influence.

+ ^h.Global experiment on model parameter sensitivity.^

    We used Sobol indices to quantify the total influence of parameters on the output: Sobol parameters are (variance-based). The Sobol parameter method (is based on variance decomposition). We first used Saltelli sampling to select a large number of parameter combinations within the parameter range as experimental samples. For each sample, the DDE model was run, the output signal was recorded, and the $S1$ (influence of a single parameter on the output), $ST$ (total influence from single parameters to two-parameter interactions up to all parameter interactions), and $S2$ (interaction strength of two parameters) for each parameter were calculated through Sobol decomposition. Based on this, we evaluated the importance of each parameter and their interactions.

The Python project was the main tool for building, debugging the model, and generating modeling results. We created the graph data structure file `graph.py` and, based on this, performed Python coding for the deterministic model, Monte Carlo method random noise simulation and LOD/LOQ calculation, parameter sensitivity analysis, and plotting, etc.

![](https://static.igem.wiki/teams/5924/drylab/ht.webp)
![](https://static.igem.wiki/teams/5924/drylab/ml.webp)

^ The interactive platform we designed has powerful and robust functionality, performance, and a user-friendly, intuitive human-computer interface. Its core functionality revolves around the dual-view of simulation results / model editor, supporting interactive editing of nodes and reaction equations, dynamic parameter configuration, `JSON` format model import/export, and still integrates noise simulation and LOD/LOQ calculation modules, capable of simulating the impact of different noise levels on the signal pathway. ^

On the html, you can easily built the model you want by changing the parameter, and quickly export your results. For example, you can change the number of two strains to optimize the ratio of them. What's more it is a fascinating platform for us to debug our program using the data provided by future experiments

[https://github.com/caoruoxi802-bot/the-plat-form/blob/main/Interactive%20Signal%20Path%20Simulation%20Platform%20(Advanced%20Noise%20Model%20Edition).html](https://github.com/caoruoxi802-bot/the-plat-form/blob/main/Interactive%20Signal%20Path%20Simulation%20Platform%20(Advanced%20Noise%20Model%20Edition).html)

Click the website to try our html!


^8.Results Analysis^

+ ^Analysis of the Deterministic DDE Model^

    To simplify the analysis, we first studied the dose-response relationship between BoNT/A concentration and signal intensity when the ratio of upstream to downstream yeast numbers was 1:1.

    Over a concentration range from very small to very large (even exceeding the biological system's limit), the BoNT/A concentration and signal intensity overall exhibit a Hill equation relationship, which is similar to reports in multiple similar references, and the Hill coefficients are relatively close (n = 1.1 [1], n = 1 (this model)). Considering the characteristics of α-factor binding to the Ste2 receptor (Hill\ n = 1), this suggests that on a large scale, the signal intensity is mainly determined by the binding of α-factor and Ste2 receptors, while the downstream signaling pathway, metabolic pathway, and upstream BoNT/A cleavage do not show significant Hill effects, approximately presenting a linear relationship, which is also consistent with our model's assumptions.
    
![](https://static.igem.wiki/teams/5924/drylab/bonta-dose-response-curve-28800s-ration-1-18to9.webp)



    In the concentration range more likely for actual biological systems (e.g., [10^{-18} M, 10^{-5} M]), the BoNT/A–signal curve shows a trend of first rising, then falling, approaching a plateau, and then continuing to rise. Preliminary analysis and estimation suggest that the declining trend and plateau characteristics may both be caused by the mutual cancellation between the increase in betaxanthin production with α-factor concentration and the increase in intermediate product degradation rate with α-factor concentration.

![](https://static.igem.wiki/teams/5924/drylab/bonta-dose-response-curve-28800-ration-1-18to-3.webp)

    Studying the working range of [10^{-18} M, 10^{-10} M] (the exact working range needs to be determined after the code results are run), it was found that the relationship between BoNT/A concentration and signal intensity again shows approximately parabolic characteristics. This may indicate that in the working range, the effects of degradation and other factors on the sensor's performance must be considered, providing important clues for the next step of research.

![](https://static.igem.wiki/teams/5924/drylab/bonta-dose-response-28800s-ration-1-18-8.webp)

+ ^Optimization of Initial Upstream and Downstream Yeast Numbers for Experiments^

![](https://static.igem.wiki/teams/5924/drylab/signal-response-curve-at-different-ration-cell-8h-all-timepin-angle1.webp)

    In the scan of the initial ratio of upstream to downstream yeast for a specified 8-hour detection time, we can easily summarize that:
    
    - When the ratio of the two strains is relatively large (e.g., greater than 10^{-3}), 8 hours will be insufficient for the betaxanthin signal to reach its current maximum even at sufficiently high BoNT/A concentrations. This may indicate an excessive number of upstream fusion proteins, where the released α-factor is not fully utilized by the downstream strain. In this concentration range, the α-factor cannot be completely captured and bound by Ste2 receptors followed by internalization and degradation within 8 hours.
    
    - When the ratio of the two strains is very small (e.g., less than 10^{-4}), the betaxanthin signal can ultimately be captured by high concentrations of BoNT/A within 8 hours, suggesting that around 10^{-4} is the threshold for full utilization of BoNT/A. Simultaneously, as the ratio of the two strains decreases, the α-factor can be fully cleaved, bound, and degraded within 8 hours even at lower BoNT/A concentrations; meanwhile, the rising limb of the curve shifts left, effectively increasing the signal intensity at lower concentrations, which will aid detection and quantification by the naked eye and instruments at low concentrations.
    
    - From the overall trend, reducing the ratio of upstream to downstream yeast will help increase the maximum signal intensity and the signal intensity at low concentrations. We speculate that the LOD and LOQ will also decrease accordingly. Since a ratio below 10^{-5} is meaningless for practical experiments, the optimal number ratio of upstream to downstream yeast can be set to 10^{-5}.

+ ^Single Parameter Sensitivity Analysis^

    Simulation results indicate that parameter sensitivity itself is related to the evaluation time and the set variation rate. It is certain that due to the chain of reactions and the set delay effects, the parameter sensitivity coefficients also have temporal stages. The variation rate mainly affects the absolute value of the parameter sensitivity to some extent, with little impact on the overall trend.

    Preliminary analysis is also relatively obvious: in the early stages of testing, the activation level of downstream reactions is small, and the signal intensity is mainly determined by upstream delay times and parameters; as the reaction time extends to 8 hours, downstream intermediate products increase significantly, and ultimately the importance of downstream parameters gradually increases; simultaneously, because the concentrations of most substances increase overall, the stability of the signal system itself also improves systematically, and the absolute values of various parameters gradually converge from 2 hours to 8 hours.

![](https://static.igem.wiki/teams/5924/drylab/single-parameter-sensitivity-analysis.webp)

    We analyzed the parameter sensitivity values at 8 hours for various variation rates and found:
    
    - From the ranking of the absolute values of sensitivity, we can Screen out the key parameters ($EC_{50}$, $k_{transcription-pGAL1-max}$, $kcat_{MjDOD}$, etc.) and corresponding processes (α-factor binding to Ste2 receptor, pGAL1 promoter transcription maximum, MjDOD catalytic turnover, etc.); in contrast, some processes are also significantly unimportant (such as protein and betaxanthin degradation, time delays, etc.).
    
    - The change in parameter sensitivity with different variation rates also reflects the importance of parameters and processes. We observed that at variation rates of 10\%, 30\%, and 45\%, parameters like $EC_{50}$, k_{transcription-pGAL1-max}, and k_{betalamicacid-to-betaxanthins} ranked at the top of parameter sensitivity, further indicating their importance. Meanwhile, as the parameter variation rate increases, the importance of Km_bonta continuously rises, suggesting the importance of confirming the kinetic parameters of BoNT/A under conditions of high parameter uncertainty.
    
    - At 2.0h, 4.0h, and 8.0h under different parameter variation rates, the types of positively and negatively correlated parameters are basically constant. However, the α-factor degradation half-life $half\_life\_alphafactor $ shows a phenomenon of negative correlation in the short term and positive correlation in the long term, but the specific cause of this change needs further exploration and discussion. Currently, there is no good conclusion, and this will be our next task.

Looking across the various set evaluation times and parameter variation rates, some parameters show obvious high sensitivity, and they are the focus of our next step in parameter investigation and experimental verification.

+ ^Global Parameter Sensitivity Analysis^

    Although single-parameter sensitivity analysis is relatively simple, it is difficult to discover the impact of interactions between parameters on experimental results. Moreover, the "single-parameter sensitivity analysis" algorithm explored above is relatively simple and cannot fully reflect the impact of parameter changes on the results. We performed a global parameter sensitivity analysis using the Sobol algorithm, focusing on capturing the impact of interactions between parameters on experimental results.

    Because the Sobol algorithm itself is based on statistical variance decomposition, its mathematical essence is different from the linear analysis of the aforementioned single-variable sensitivity analysis. Therefore, its results under the same conditions may also differ from the previous method, but this also provides us with another perspective for analysis.

    We analyzed the global parameter sensitivity analysis results at 8 hours with parameter variation rates Δ of $20\%$ and $30\%$, respectively, as shown in the figures below. We found that the two situations are relatively similar, so we merged the analysis and found many conclusions similar to the single-parameter analysis:

    It must be specially declared that the following analysis no longer includes time delays and the reaction coefficient of the spontaneous betaxanthin reaction in the scope of parameter sensitivity research. Due to time constraints, only the following few parameter sensitivities were studied, which differs from the above.

![](https://static.igem.wiki/teams/5924/drylab/global-parameter-sensitivity-analysis-8h-var0-30.webp)

![](https://static.igem.wiki/teams/5924/drylab/heatmap-parameter-interaction-8h-var0-30.webp)

![](https://static.igem.wiki/teams/5924/drylab/parameter-sensitivity-contribution-analysis-8h-var0-30.webp)

    - Under long time (8h) and high variation rate (30\%), the Km value of BoNT/A becomes one of the most sensitive parameters, but its sensitivity is similar to most nodes, mainly manifested through interactions with other nodes.
    
    - Important parameters such as ration_cell, $k_{transcription-GAL1-max$, $EC_{50}$ still rank the forefront.
    
    But at the same time, we note that by considering interactions, we observe that the interactions between parameters reflect more meaningful properties:
    
    - $Km_{bonta}$ and $kcat_{bonta}$ are both important kinetic parameters of BoNT/A. Their close coupling in parameter sensitivity naturally reflects the connection of enzyme kinetic properties.
    
    - We understand that in complex and highly coupled biochemical systems, ST is often several times higher than S1, which is similar to our analysis results, reminding us to pay greater attention in the next step.
    
    - Surprisingly, $Km_{bonta}$ had the strongest interaction with k_{translation-mjdod}$ in terms of parameter sensitivity. There is currently no good explanation for this, awaiting further analysis and research.

^9. Guidance to Wet Lab^

The above model results and findings have already had an important impact on experimental work:

- We determined that the optimal input ratio of the two yeast strains (upstream and downstream) before the reaction, which minimizes the Limit of Detection (LOD) at an 8-hour reaction time, can be as low as 0.0001. Lower ratios have no significant meaning for experiments, thus determining the best initial reaction conditions.

- The model's prediction of the relationship between signal, time, and BoNT/A concentration provides clues for controlling reaction time and narrows down the possible concentration range of initial substrates in real experiments and simulation experiments. Our next simulation experiments can focus on the region near the predicted LOD and LOQ, avoiding scanning over an extremely large concentration range.

- The model's analysis of local and global parameter sensitivity implies that it will be necessary to strengthen experiments and data collection on sensitive parameters in the next step, guiding the focus of future experimental work.

We also shared the modeling process and results synchronously with the entire team through intuitive forms like `html` web pages, promoting the exchange and dissemination of knowledge within the team.

![](https://static.igem.wiki/teams/5924/drylab/ht.webp)
![](https://static.igem.wiki/teams/5924/drylab/ml.webp)

[https://github.com/caoruoxi802-bot/the-plat-form/blob/main/Interactive%20Signal%20Path%20Simulation%20Platform%20(Advanced%20Noise%20Model%20Edition).html](https://github.com/caoruoxi802-bot/the-plat-form/blob/main/Interactive%20Signal%20Path%20Simulation%20Platform%20(Advanced%20Noise%20Model%20Edition).html)

Click the website to try our platform!

Of course, experimental work itself is the source of modeling. The continuous collection and organization of experimental data and experimental intelligence have helped us determine and confirm some processes, initial conditions, and parameters:

a. We determined the initial concentration of L-DOPA; and through discussion with the Wet lab, determined the initial value of the system's amino acid concentration.

b. Through analysis of wet lab raw data, we plotted the overall absorption spectrum of the pathway from MjDOD catalyzing L-DOPA to betalamic acid and downstream, confirming that the overall maximum absorption peak is around 480nm.

c. We collaborated with the Wet Lab to further design improvement plans for future work based on the modeling results, such as adding verification experiments and improving the expression system.

The above shows the close connection between the model and experiments.

^10.Discussion& Future Work^

This model successfully derived the BoNT/A–signal intensity dose-response curves at different specific times and successfully predicted the LOD and LOQ by introducing randomness. This is the final realization and effect test of the model design. It indicates the effectiveness and shortcomings of the current design and guides the next round of optimized design. The analysis of parameter sensitivity reminds us to focus on confirming and optimizing some parameters.

Although the form of this model is similar to conventional ODE modeling of signaling and metabolic pathways, we still strived for innovation in the model: we imaginatively linked the pathway to a "graph," preserving the structural characteristics of the network while enabling quantitative simulation of dynamic processes, providing a new perspective for time-varying analysis of complex biological networks.

We introduced Monte Carlo simulation based on Gaussian random noise in key processes of the biological system. Compared with general noise simulation methods, the idea and implementation are simpler, and it is easy to flexibly adjust the noise level. Layered noise processing also conforms to the laws of biological systems. The combination of signal pathways and graph theory introduced Dijkstra's algorithm into signal pathway research and also gave rise to our application of Delay Differential Equations (DDEs).

We developed a web simulation platform with good scalability and intuitiveness, which helps promote knowledge dissemination. Its degree of personalization can also be higher than conventional software like COPASI, while reducing the learning cost for beginners to get started with pathway modeling.

Of course, no modeling of biological systems can be called completely correct, and our model is no exception. Due to the great uncertainty of the parameters themselves, the simplification of linear equation expressions, the neglect of secondary processes, and the inherent chaos of biological systems, our model must have many defects, such as insufficient capture of important processes in the pathway and large discrepancies in predicted data.

Furthermore, due to the rush in research time, our model also has many unfinished tasks. Our research on the optimal ratio of the two strains and parameter sensitivity is still in its early stages, and we have not fully explored the connotations and deep relationships behind them. Although generally consistent, our linear model also showed some inconsistencies with reality in simulations. We will bravely face challenges and difficulties, constantly evolving with the help of the "Design-Build-Test-Learn" research flywheel.

Nevertheless, our model has already provided us with much useful information. The exploration of the optimal ratio of the two strains, the study of LOD and LOQ, and the analysis of local and global sensitivity reminds us that the subsequent work plan should focus on the above issues, such as model parameter optimization, process tuning, and initial condition optimization.

We will deeply communicate with the Wet Lab to determine the best parameter values, and by comparing with experimental results, uncovering overlooked channels, better practicing the "Design-Build-Test-Learn" cycle. We will introduce more reasonable randomness and uncertainty analysis to enhance the robustness of LOD and LOQ analysis. Driven by the dual wheels of experiment and modeling, we will refine and optimize the reliability of the project design, striving to design a stable and sensitive biosensing system. We will also try to introduce more graph theory methods into biological system modeling, providing unique research perspectives.

We will also further clarify the development direction and content of the web platform, try to integrate more powerful functions, and inherit the intuitive and easy-to-understand nature, turning it into a more inclusive and practical learning platform, making it easy to get started  biological system modeling.

^11. References^

[1] V. Iyer, & K. Struhl, Absolute mRNA levels and transcriptional initiation rates in Saccharomyces cerevisiae., Proc. Natl. Acad. Sci. U.S.A. 93 (11) 5208-5212, https://doi.org/10.1073/pnas.93.11.5208 (1996).

[2]  Yi, H. Kitano, & M.I. Simon, A quantitative characterization of the yeast heterotrimeric G protein cycle, Proc. Natl. Acad. Sci. U.S.A. 100 (19) 10764-10769, https://doi.org/10.1073/pnas.1834247100 (2003).

[3] Sherman F. (2002). Getting started with yeast. Methods in enzymology, 350, 3–41.

[4] [PIR1 | SGD](https://www.yeastgenome.org/locus/S000001647)

[5] Fan, C., He, N., & Yuan, J. (2024). Cascaded amplifying circuit enables sensitive detection of fungal pathogens. Biosensors and Bioelectronics, 250, 116058, https://doi.org/10.1016/j.bios.2024.116058

[6] Chou YC, Shih CI, Chiang CC, Hsu CH, Yeh YC. Reagent-free DOPA-dioxygenase colorimetric biosensor for selective detection of L-DOPA. Sensors and Actuators B: Chemical. 2019 Oct 15;297:126717. https://doi.org/10.1016/j.snb.2019.126717

[7] Mizanur RM, Stafford RG, Ahmed SA (2014) Cleavage of SNAP25 and Its Shorter Versions by the Protease Domain of Serotype A Botulinum Neurotoxin. PLoS ONE 9(4): e95188. https://doi.org/10.1371/journal.pone.0095188

[8] calculated by protparam https://web.expasy.org/cgi-bin/protparam/protparam

[9] https://www.tipbiosystems.com/wp-content/uploads/2023/12/AN102-Yeast-Cell-Count_2019_03_17-1.pdf

[10] Hsu, C., Scherrer, S., Buetti-Dinh, A. et al. Stochastic signalling rewires the interaction map of a multiple feedback network during yeast evolution. Nat Commun 3, 682 (2012). https://doi.org/10.1038/ncomms1687

[11] [NanoDrop™ One/OneC Microvolume UV-Vis Spectrophotometer ND OneC Spec, no Wi-Fi | Request for Quote](https://www.thermofisher.cn/order/catalog/product/cn/zh/840-317500)

[12] Girod, PA., Zryd, JP. Secondary metabolism in cultured red beet (Beta vulgaris L.) cells: Differential regulation of betaxanthin and betacyanin biosynthesis. Plant Cell Tiss Organ Cult 25, 1–12 (1991). https://doi.org/10.1007/BF00033905

[13] Gao, X., Lanz, M., Grosely, R., Cremer, J., Puglisi, J., & Skotheim, J. M. (2025). Yeast growth is controlled by the proportional scaling of mRNA and ribosome concentrations. *arXiv preprint arXiv:2508.14997.*, https://arxiv.org/abs/2508.14997

[14] [DOD - 4,5-DOPA dioxygenase extradiol - Mirabilis jalapa (Garden four-o'clock) | UniProtKB | UniProt](https://www.uniprot.org/uniprotkb/B6F0W8/entry)

[15] Y. Arava, Y. Wang, J.D. Storey, C.L. Liu, P.O. Brown, & D. Herschlag, Genome-wide analysis of mRNA translation profiles in Saccharomyces cerevisiae, Proc. Natl. Acad. Sci. U.S.A. 100 (7) 3889-3894, https://doi.org/10.1073/pnas.0635171100 (2003).

[16] (Same as [14])

[16] Contreras-Llano LE, Guerrero-Rubio MA, Lozada-Ramírez JD, García-Carmona F, Gandía-Herrero F. 2019. First Betalain-Producing Bacteria Break the Exclusive Presence of the Pigments in the Plant Kingdom. mBio 10:10.1128/mbio.00345-19. https://doi.org/10.1128/mbio.00345-19

[17] Chan, L. Y., Mugler, C. F., Heinrich, S., Vallotton, P., & Weis, K. (2018). Non-invasive measurement of mRNA decay reveals translation initiation as the major determinant of mRNA stability. Elife, 7, e32536, https://doi.org/10.7554/eLife.32536

[18] Aymoz, D., Solé, C., Pierre, J. J., Schmitt, M., de Nadal, E., Posas, F., & Pelet, S. (2018). Timing of gene expression in a cell‐fate decision system. Molecular systems biology, 14(4), e8024.

[19] Steinbrecht, D., Minia, I., Milek, M., Meisig, J., Blüthgen, N., & Landthaler, M. (2024). Subcellular mRNA kinetic modeling reveals nuclear retention as rate-limiting. Molecular Systems Biology, 20(12), 1346-1371, https://doi.org/10.1038/s44320-024-00073-2

[20] Cruz, S., Checa, N., Tovar, H., Cejudo-Bastante, M. J., Heredia, F. J., & Hurtado, N. (2024). Semisynthesis of Betaxanthins from Purified Betacyanin of Opuntia dillenii sp.: Color Stability and Antiradical Capacity. Molecules, 29(9), 2116. https://doi.org/10.3390/molecules29092116

[21]Lamphier, M. S., & Ptashne, M. (1992). Multiple mechanisms mediate glucose repression of the yeast GAL1 gene. *Proceedings of the National Academy of Sciences*, *89*(13), 5922-5926.https://doi.org/10.1073/pnas.89.13.5922

# Optimize the ratio of Substrate-acting Strain and Signal-transferring Strain

^1. Introduction^

The co-culture of Substrate-acting Strain and Signal-transferring Strain is the foundation of our system's communication mechanism. Substrate-acting strain secrete α-factor, while Signal-transferring strain detect and respond to it through the Ste2p receptor pathway. To ensure effective communication, the ratio between Substrate-acting Strain and Signal-transferring Strain is crucial: too few Substrate-acting Strains lead to insufficient signaling, while too many Substrate-acting Strains cause signal saturation and resource waste.

Therefore, we constructed ^three^ mathematical models to explore how the Substrate-acting Strain-to-Signal-transferring Strain ratio influences α-factor distribution and receptor activation.^Since the quantities of the two strains in the initial system have reached to maximum , we do not consider the competition during yeast reproduction in the model.^ By quantitatively describing the balance between α-factor production, diffusion, and uptake, our model aims to identify the optimal ratio that maximizes signaling efficiency. This provides both a theoretical explanation for the experimental observations and guidance for adjusting proportions in the wet lab.

# Model 1: 

^1. Model Assumptions^

    a. The spatial transport of α-factor between Substrate-acting Strain and Signal-transferring Strain can be approximated as one-dimensional diffusion, since the main direction of α-factor transfer is from Substrate-acting Strain cells to signal cells across the medium.

    b. α-factor degradation within the bulk medium is negligible on the time scale of interest; removal occurs primarily via boundary uptake.

    c. The system operates under steady temperature and homogeneous medium conditions, so the diffusion coefficient $D$ is constant.

^2. Equation Construction^

The α-factor concentration $C(x,t)$ is modeled by a one-dimensional diffusion equation over the domain $x\in[0,L]$:\[\frac{\partial{C}}{\partial{x}}=D\frac{\partial^2{C}}{\partial{x^2}}\]
The Substrate-acting Strain at $x=0$ secrete α-factor as a flux source $R$,while the Signal-transferring Strain at $x=L$ absorb it through a constant uptake flux $J_0$. Thus, the boundary conditions are:\[-D\frac{\partial{C}}{\partial{x}}(0,t)=R,\quad -D\frac{\partial{C}}{\partial{x}}(L,t)=R\]
The production rate follows a Michaelis–Menten form,\[R=\frac{K_{cat}[BoNT][SNAP]_{total}}{K_m+[SNAP]_{total}}\]while the uptake rate is approximated as a constant $J_0=k_3N_{recv}$,independent of local concentration.

^3. Parameter List^

| Parameter | Meaning | Value | Unit | Reference |
|:---:| :---:| :---:| :---:| :---: |
|$D$|Diffusion coefficient of α-factor in the medium|$2.18\times10^{-10}$|$$m^2/s$$| [3] |
|$C(x,t)$|Concentration of α-factor at position $x$ and time $t$|$mol/m^3$| |
|$K_{cat}$|Turnover number of BoNT acting on SNAP-25|$60$|$s^{-1}$| [1] |
|$K_m$|Michaelis constant of SNAP-25 cleavage by BoNT|$1.62\pm0.6$|$\mu M$| [1] |
|$k_3$|Effective constant uptake coefficient at Signal-transferring Strain|$1\times 10^{-7}$|$mol/(cell\cdot s)$| estimated |


^4. Equation Solve^

This PDE is solved numerically using the finite element method (FEM). The weak form yields a semi-discrete system\[\textbf{M}\frac{d\textbf{u}}{dt}+D\textbf{Ku}=\textbf{F}(t)-\textbf{G}(t),\]where $\textbf{M}$ and $\textbf{K}$ are the standard FEM mass and stiffness matrices, and $\textbf{F}$,$\textbf{G}$ represent the production and uptake fluxes at the boundaries. Time integration is performed using the implicit Euler scheme:\[(\textbf{M}+\Delta tD\textbf{K})\textbf{u}^{n+1}=\textbf{Mu}^n+\Delta t(\textbf{F}^{n+1}-\textbf{G}^{n+1})\]The cumulative absorbed amount is then obtained by integrating $J_0$ over time.

^5. Result^

We scanned the Substrate-acting Strain fraction $p$ from 0.01 to 0.99 and computed the cumulative α-factor uptake $J_2$ by the Signal-transferring Strain using our 1D diffusion model. The results show that the optimal ratio requires a significantly higher number of Signal-transferring Strain compared to Substrate-acting Strain. In other words, effective communication is achieved when Signal-transferring Strain dominate the population, allowing them to efficiently capture α-factor produced by fewer Substrate-acting Strains.

![](https://static.igem.wiki/teams/5924/drylab/ratio.webp)

^6. Conclusion^

The model suggests that Signal-transferring Strain should outnumber Substrate-acting Strain to maximize cumulative uptake, with Substrate-acting Strains forming only a small fraction of the population. This ensures sufficient α-factor signaling without saturation. While some parameters, like the uptake coefficient $k_3$ , are estimated, the qualitative trend is robust: increasing Signal-transferring Strain improves overall signaling efficiency. 

# Model 2:

To verify the correctness of the model, we also used the modeling program employed in the entire process simulation of the system to model the proportioning.

^1. Method^

The total number of the two types of strains remains unchanged.(2e11 cells/L,[5])The concentration of Signal-transferring strain decides the concentration of ste2 receptor,while the concentration of Substrate-acting strain decides the  concentration of Snap25.

More Signal-transferring strain leads to higher possibility of ste2 receptors coming into contact with alpha factors ,thus more more betaxanthin would be produced  which means the color are easier to recognize. 

Meanwhile more Substrate-acting strain leads to more snap 25, thus more α-factor would be freed which means more ste2 would be actived.

\[C_{Signal-transferring}=C_{snap25}/x\]

In the entire process simulation of the system, the concentration of snap25 is related to that of Substrate-acting as mentioned above. Therefore, I can adjust the number of cells by regulating the concentration of snap25.

\[C_{Substrate-acting}+C_{Signal-transferring}=K_{cell}\]

\[ratio=\frac{C_{Signal-transferring}}{K_{cell}}\]

| Parameter | Meaning | Value | Unit | Reference |
|:---:| :---:| :---:| :---:| :---: |
|$C_{snap25}$|The concentration of Substrate-acting strain within the system|-|$M$| - |
|$x$|the number of snap25 on each Substrate-acting strain |3000|1| estimated |
|$C_{Substrate-acting}$|The cell concentration of Substrate-acting strain within the system|-|$$cells/L$$| - |
|$C_{Signal-transferring}$|The cell concentration of Signal-transferring strain within the system|-| $$cells/L$$|-|
|$K_{cell}$|Environmental carrying capacity for two strains in YPD|$2e11$|$cells/L$| [7] |

^2. Result^

![](https://static.igem.wiki/teams/5924/drylab/rat2.webp)

By altering the  the proportion of signal-transferring strains  in the entire strain in its system, we obtained the maximum signal intensity corresponding to a response of 50,000 seconds for different botulinum toxin concentrations. We could find that the more the signal-transferring strains ,the stronger the signal. The plot shows the same trend as the results obtained from the previous modeling, verifying the correctness of the proportioning modeling.

![](https://static.igem.wiki/teams/5924/drylab/qujian.webp)
At the same time, we also found that the the effective concentration range of botulinum toxin in the test kit vary with different ratios.We can find that the higher the proportion, the better the system performs at low concentrations (1e-16M), and the lower the proportion, the better the system performs at higher concentrations (1e-15M).



# Model 3:
In this model, we constructed a kinetic network for α-factor release and binding, incorporating competition between the two strains. In the signal-transferring strain, the native Gal4 promoter in the genome was deleted. According to previous studies, deletion of the genomic Gal4 promoter markedly reduces yeast growth rate [6]. Therefore, the substrate-acting strain exhibits a faster growth rate than the signal-transferring strain. 

^1.Method^

The dynamics of the system are captured by the following set of ordinary differential equations (ODEs):

\[\frac{d[U]}{dt} = r_U \cdot [U]\cdot (1-\frac{[U]+[D]}{C_{tot}})\]

\[\frac{d[D]}{dt} = r_D \cdot [D]\cdot (1-\frac{[U]+[D]}{C_{tot}})\]

\[\frac{d[A]}{dt} = v_{cut}- k_{on} [A][R_f]+k_{off}[R_b]\]

\[\frac{d[R_b]}{dt} = k_{on}[A][R_f]- k_{off} [R_b]-k_{int}[R_b]\]

\[\frac{d[R_f]}{dt} = -k_{on}[A][R_f]+ k_{off} [R_b]+R_{0c}\cdot \frac{d[D]}{dt}\]

\[\frac{d[F_u]}{dt} = -v_{cut}+F_{0c}\cdot \frac{d[U]}{dt}\]

\[v_{cut}=\frac{k_{cat}\cdot [F_u]\cdot Toxin}{K_m+[F_u]}\]

\[\frac{d[S]}{dt} = k_s[R_b]\]

^2. Variable list ^

| Symbol | Meaning | Initial Value |
| :-------- | :------------- | :-------------- |
| $[U]$ | Concentration of upstream strain (substrate-acting strain) | $f \cdot C_{tot} \cdot p$ uM |
| $[D]$ | Concentration of downstream strain (signal-transferring strain) | $(1-f)\cdot C_{tot} \cdot p$ |
| $[A]$ | Concentration of α fatcor | 0 |
| $[R_b]$ | Concentration of bound Ste2 receptors | 0 |
| $[R_f]$ | Concentration of unbound Ste2 receptors | $R0_c\cdot [D]$ |
| $[F_u]$ | Concentration of cleavable fusion protein | $F0_c\cdot [U]$ |
| $v_{cut}$ | Cleavage rate | – |
| $[S]$ | Relative signal intensity | 0 |


^3. Parameter list ^

| Symbol | Meaning | Value | Unit | References |
| :-------- |:-----------------|:-------------------------|:---------------| :-------- |
| $f$ | Proportion of substrate-acting strain in the total cell population | Variable | - | - |
| $C_{tot}$ | Maximum yeast concentration | $2\times10^8$ | cells/mL | [7] |
| $p$ | Initial fraction of the maximum yeast concentration | 0.5 | - | assumed |
| $r_U$ | Growth rate of substrate-acting strain | 0.0001 | s-1 | estimated |
| $r_D$ | Growth rate of signal-transferring strain | 0.00005 | s-1 | estimated |
| $F0_c$ | Number of fusion proteins per substrate-acting strain cell | $1.58 \times 10^{-12}$ | nmol/cell | estimated from gene expression model|
| $R0_c$ | Number of Ste2 receptors per signal-transferring strain cell | $1.33 \times10^{-11}$ | nmol/cell | [8] |
| $k_{cat}$ | Catalytic rate of botulinum toxin cleavage | 76 | s-1 | [9] |
| $K_m$ | Michaelis constant of botulinum toxin | 33.34 | uM | [9] |
| $Toxin$ | Concentration of botulinum toxin | $1\times10^{-6}$ | uM | |
| $k_{on}$ | Association rate of Ste2 and α-factor | $1\times10^{-3}$ | uM-1s-1 | assumed |
| $K_d$ | Dissociation constant of Ste2 and α-factor | $5.6\times10^{-3}$ | uM | [10] |
| $k_{off}$ | Dissociation rate of Ste2–α-factor | $k_{on} \cdot K_d$ | s-1 | assumed |
| $k_{int}$ | Receptor internalization rate | $1\times10^{-3}$ | s-1 | assumed |
| $k_s$ | Signal transduction rate | $1\times10^{-3}$ | s-1 | assumed |
| $t_{end}$ | Simulation end time | 6000 | s | - |
| $dt$ | Integration time step | 0.5 | s | - |

^4. Result ^

The results show that the signal intensity reaches its maximum when (f = 0.5), corresponding to a 1:1 ratio of substrate-acting strain to signal-transferring strain.
![](https://static.igem.wiki/teams/5924/drylab/ode-ratio.webp)


# Discussion:

During the modeling of signal transmission efficiency in a two-strain co-culture system, we employed three distinct approaches to simulate the binding of α-factor to the Ste2 receptor and downstream signal activation.

In the first model, we assumed each signal-transferring cell produces a fixed signal, independent of the surrounding α-factor concentration: 

\[
J_0 = k_3 N_{recv}
\]

Under this assumption, the total signal depends only on the product of single-cell signal and signal-transferring strain number. As a result, the model fails to accurately reflect the contribution of the substrate-acting strain, particularly when the free α-factor concentration is low (i.e., at low substrate-acting strain proportions).

Next, the second model incorporates concentration-dependent signaling, using an EC50-like formulation to describe signal strength:

\[
\theta = \frac{c_{\alpha\text{-factor-free}}}{EC50_{\text{ste2}} + c_{\alpha\text{-factor-free}}}
\]

This model reflects signal enhancement with rising α-factor levels and includes a saturation limit via EC50. Since θ is restricted between 0 and 1, the total signal (i.e., the product of θ and signal-transferring cell number) remains strongly dependent only on the number of signal-transferring cells.

Finally, the third model simulates signal activation using ligand–receptor binding kinetics:

\[
\frac{d[R_b]}{dt} = k_{\text{on}}[A][R_f] - k_{\text{off}} [R_b] - k_{\text{int}} [R_b]
\]

In this model, if the total cell number is $N$ and the proportion of producer cells is $p$, the amount of α-factor produced is approximately proportional to $pN$. The uptake rate is roughly proportional to both the α-factor concentration and the number of signal-transferring cells. As a result, regardless of other parameter changes, the total uptake follows a quadratic function proportional to $p(1-p)$, reaching its maximum near (p = 0.5), indicating a balanced contribution from both upstream and downstream strains.

# Advice:

Based on our modeling results, we recommend that future experiments explore two approaches: (1) maintaining a high proportion of signal-transferring cells relative to substrate-acting cells, and (2) keeping the ratio of signal-transferring to substrate-acting cells close to 1:1. The optimal strategy should be further validated through wet-lab experiments.



# References

[1]Unique Substrate Recognition by Botulinum Neurotoxins Serotypes A and E *Chen, Sheng et al. Journal of Biological Chemistry, Volume 281, Issue 16, 10906 - 10911
[2]Moore TI, Chou C-S, Nie Q, Jeon NL, Yi T-M (2008) Robust Spatial Sensing of Mating Pheromone Gradients by Yeast Cells. PLoS ONE 3(12): e3865.
[3]Jahn M, Mölle A, Rödel G, Ostermann K. Temporal and spatial properties of a multi-cellular amplification system based on signal molecule diffusion. Sensors (Basel). 2013 Oct 25;13(11):14511-22.
[4]Chen W, Nie Q, Yi TM, Chou CS. Modelling of Yeast Mating Reveals Robustness Strategies for Cell-Cell Interactions. PLoS Comput Biol. 2016 Jul 12;12(7):e1004988.
[5]Sherman F. Getting started with. Methods Enzymol. 2002;350:3-41. 
[6] Giaever, G., Chu, A., Ni, L. et al. Functional profiling of the Saccharomyces cerevisiae genome. Nature 418, 387–391 (2002). 
[7] Sherman F. (2002). Getting started with yeast. Methods in enzymology, 350, 3–41.
[8] Jenness, D. D., Burkholder, A. C., & Hartwell, L. H. (1986). Binding of alpha-factor pheromone to Saccharomyces cerevisiae a cells: dissociation constant and number of binding sites. Molecular and cellular biology, 6(1), 318–320.
[9] Mizanur, R. M., Stafford, R. G., & Ahmed, S. A. (2014). Cleavage of SNAP25 and its shorter versions by the protease domain of serotype A botulinum neurotoxin. PloS one, 9(4), e95188.
[10] Bush, A., Vasen, G., Constantinou, A., Dunayevich, P., Patop, I. L., Blaustein, M., & Colman-Lerner, A. (2016). Yeast GPCR signaling reflects the fraction of occupied receptors, not the number. Molecular systems biology, 12(12), 898.
