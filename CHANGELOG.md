12.08.2016:
* Added a CHANGELOG

13.05.2019:
* first update in CHANGELOG
* only documentation of principal version changes in CHANGELOG

LeMonADE 1.0
* bare copy of LeMonADE library from ITP (Leibniz-Institute of Polymer Research Dresden IPF) with basic functionalities

LeMonADE 2.0
* adding abstract system setup Updater: UpdaterAbstractCreate
* adding the nearest neighbor interaction Features: FeatureNNInteractionSc and FeatureNNInteractionBcc

LeMonADE 2.1
* adding external potential Feature: FeatureSpringPotentialTwoGroups
* adding dynamic bond changing: FeatureConnectionSc
* modify FeatureAttribute: now uses a template parameter `FeatureAttribute<>`, default=`int32_t`
* improvement of utilities: DistanceCalculation
* introduce continuous integration: CircleCi

LeMonADE 2.2
* add moveable Tags
* add an updater for a swelling box
* add analyzers for the calculation of the mean-square displacement
* add various features and moves to enable reactive and reversible bonds

LeMonADE 2.2.2
* HotFix for default R250 values
* ColdFix for implementation of RNG Mersenne Twister 
* ColdFix for ResultFormattingTools

LeMonADE 2.3.0
* add *MonomerInteractionTag* for feature FeatureNNInteractionSc
* add additional lattice for interactions in FeatureNNInteractionSc