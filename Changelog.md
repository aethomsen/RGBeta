### v.1.1.3 (2023-01-30)
- Fixed compatibility issue with Mathematica 13.2.

### v.1.1.2 (2022-05-19)
- Implemented mixed terms in AnomalousDimension and AnomalousDimTerm to capture the mixing between multiple fields of similar quantum numbers. Simply replace field with {field1, field2} in such cases.
- Similar for UpsilonTerm and UpsilonFunction.

### v.1.1.1 (2022-01-17)
- Fixed compatibility issue with Mathematica 13.0.

## v.1.1.0 (2021-10-12)
- Added 4-loop gauge beta function.
- Added 3-loop Yukawa beta function.
	- New option *FlavorImproved->True* fixes the use of the flavor-improved beta function (True by default).   
- Added UpsilonTerm and UpsilonFunction to return upsilon of a field (in the manner of AnomalousDimTerm/AnomalousDimension).

### v.1.0.3 (2021-06-16)
- Fixed bug with antisymmetric Yukawa couplings making the projection operator vanish.
- Added options SymmetricIndices and AntisymmetricIndices to specify coupling symmetries and simplify associated matrix contractions.

### v.1.0.2 (2021-04-21)
- Fixed bug with conjugation of generators in GroupInvariants.
- Implemented ordering of del[rep, a ,b] for performance.

### v.1.0.1 (2021-02-11)
- Implemented a warning in BetaFuntion and BetaTerm for quartic couplings when their projectors mix with other quartic couplings.
