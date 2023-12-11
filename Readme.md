# Arts2024ContractionModel

<https://github.com/CircAdapt/Arts2024ContractionModel>

Source code of the CircAdapt model of heart and circulation, in which is the mechanochemical-based contraction model is implemented, as described in the manuscript: “Translating myosin-binding protein C and titin abnormalities to whole-heart function using a novel calcium-contraction coupling model” by Theo Arts1, PhD, Aurore Lyon1, PhD, Tammo Delhaas1, PhD, Diederik WD Kuster2, PhD, Jolanda van der Velden2, PhD, Joost Lumens1, PhD.
1 Department of Biomedical Engineering, Cardiovascular Research Center Maastricht (CARIM), Maastricht University, 6200MD Maastricht, The Netherlands.
2 Department of Physiology, Amsterdam University Medical Center, 1081HZ Amsterdam, The Netherlands.


# Code of the implementation of the new contraction model in CircAdapt

The code is written in MATLAB (MathWorks, version R2023a). The model of calcium-contraction coupling is implemented by the files 'SarcSf.m' and 'SarcDot.m'. In 'SarcSf.m', the sarcomere is represented by a sarcomere length and a series elastic element. These data are used to calculate mechanical equilibria in the CircAdapt model of heart and circulation. In 'SarcDot.m' the time derivatives of the state variables unloaded sarcomere length P.Patch.Lsi (=Ls in Eq.1 of the manuscript) and amount of attached cross-bridges P.Patch.Xb (=SXbTot in Eq.10 of the manuscript).

To run the CircAdapt model of  the heart and circulation, open MATLAB in folder CircAdaptXC07. Execute the script 'CircAdaptMain', selecting the reference state by pressing [R].

In the manuscript, parameters bc and ? are represented in the MATLAB-code by P.Patch.LnbC= ln(bc) and P.Patch.gTit= ?, respectively. Parameter bc indicates the factor of suppression for Ca sensitivity of the C-zone of the thick filament. Parameter ? indicates the effect of titin lengthening on the availability of cross-bridge heads to attach to the thin filament.
Time dependent stress and strain are represented by P.Patch.Sf and P.Patch.Ef, respectively.
Sarcomere length and normalized equilibrium calcium concentration are represented by P.Patch.Ls (=2 hs including series elasticity) and P.Patch.Ca (= Caeq in the manuscript).

The sequence of 5 parameters in P.Patch follow the sequence: left atrium, right atrium, left ventricular free wall, septum and right ventricular free wall.


# Abstract of the related manuscript

Mutations in cardiac myosin-binding protein C (cMyBP-C) or titin may respectively lead to hypertrophic (HCM) or dilated (DCM) cardiomyopathies. The mechanisms leading to these phenotypes remain unclear because of the challenge of translating cellular abnormalities to whole-heart and system function.
We developed and validated a novel computer model of calcium-contraction coupling incorporating the role of cMyBP-C and titin based on the key assumptions: 1) tension in the thick filament promotes cross-bridge attachment mechanochemically, 2) with increasing titin tension, more myosin heads are unlocked for attachment, and 3) cMyBP-C suppresses cross-bridge attachment.
Simulated stationary calcium-tension curves, isotonic and isometric contractions, and quick release agreed with experimental data. The model predicted that a loss of cMyBP-C function decreases the steepness of the calcium-tension curve, and that more compliant titin decreases the level of passive and active tension and its dependency on sarcomere length. Integrating this cellular model in the CircAdapt model of the human heart and circulation showed that a loss of cMyBP-C function resulted in HCM-like hemodynamics with higher left ventricular end-diastolic pressures and smaller volumes. More compliant titin led to higher diastolic pressures and ventricular dilation, suggesting DCM-like hemodynamics. 
The novel model of calcium-contraction coupling incorporates the role of cMyBP-C and titin. Its coupling to whole-heart mechanics translates changes in cellular calcium-contraction coupling to changes in cardiac pump and circulatory function and identifies potential mechanisms by which cMyBP-C and titin abnormalities may develop into HCM and DCM phenotypes. This modeling platform may help identify distinct mechanisms underlying clinical phenotypes in cardiac diseases.


