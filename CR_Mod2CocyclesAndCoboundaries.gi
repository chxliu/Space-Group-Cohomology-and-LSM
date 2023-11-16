
#####################################################################
#####################################################################

CR_Mod2CocyclesAndCoboundaries:=function(arg)
local
	R, n, toggle, Dimension, Boundary, 
	M1, M2, row, sol,
	kerdim, imgdim, cohdim, Mod2Cohomologydim,
    BasisKerd1, BasisImaged2, Rels, CobandCoc,
	#Smith, SmithRecord, TorsionCoefficients,
	ColMat, InvColMat,
	RemoveRowsMat, InsertRowsList,
    GF2ToZ,
	CycleToClass, ClassToCycle,
	i, j, k, x, sum;


R:=arg[1];
n:=arg[2];
Dimension:=R!.dimension;
Boundary:=R!.boundary;
toggle := true;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

#####################################################################
GF2ToZ:=function(v)
local v0,k;

v0:=[];
for k in [1..Length(v)] do
    if v[k] = 0*Z(2) then
    v0[k]:=0;
    else
        v0[k]:=1;
    fi;
od;

return v0;
end;
#####################################################################


	################CONSTRUCT BOUNDARY MATRICES M1 AND M2########
M1:=[];
M2:=[];

#M1 is Dim(n) x Dim(n-1);
#M2 is Dim(n+1) x Dim(n);

for i in [1..Dimension(n)] do
row:=[];
        for j in [1..Dimension(n-1)] do
        sum:=0;
                for x in Boundary(n,i) do
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:= RemInt(sum,2);
        od;
M1[i]:=row;
od;

if Dimension(n+1)>0 then
for i in [1..Dimension(n+1)] do
row:=[];
        for j in [1..Dimension(n)] do
        sum:=0;
                for x in Boundary(n+1,i) do
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:= RemInt(sum,2);
       	od;
M2[i]:=row;
od;

else

row:=[];
for j in [1..Dimension(n)] do
row[j]:=0;
od;
M2[1]:=row;
fi;
	################MATRICES M1 AND M2 CONSTRUCTED###############

#M1 is Dim(n) x Dim(n-1);
#M2 is Dim(n+1) x Dim(n);
#BasisKerd1:=LLLReducedBasis(TransposedMat(M2),"linearcomb").relations;
#BasisImaged2:=LLLReducedBasis(TransposedMat(M1)).basis;


if M2 = [ [ ] ] then
BasisKerd1:=[];
else
BasisKerd1:=BasisNullspaceModN(TransposedMat(M2),2);
fi;

if M1 = [ ] then
BasisImaged2:=[];
else
BasisImaged2:=BaseMat(TransposedMat(M1)*Z(2));
fi;

#Print(BasisKerd1);
#Print(BasisImaged2);

imgdim:=Length(BasisImaged2);
kerdim:=Length(BasisKerd1);
cohdim:=kerdim-imgdim;


CobandCoc:=[];
for i in [1..imgdim] do
    Append(CobandCoc,[BasisImaged2[i]]);
od;


if cohdim >0 then
for i in [1..kerdim] do
        if imgdim = 0 then
        Append(CobandCoc,[BasisKerd1[i]*Z(2)]);
        else
        sol:=SolutionMat(CobandCoc,BasisKerd1[i]*Z(2));
        if sol=fail then
            Append(CobandCoc,[BasisKerd1[i]*Z(2)]);
        fi;
        fi;
od;
fi;

for i in [1..kerdim] do
    CobandCoc[i]:=GF2ToZ(CobandCoc[i]);
od;


if toggle=false then 
return rec(
		#cocyclesBasis:=BasisKerd1,
		#boundariesCoefficients:=Rels,
		#torsionCoefficients:=fail,
        cocyclesBasis:=CobandCoc,
        Mod2Cohomologydim:=cohdim,
		cocycleToClass:=fail,
		classToCocycle:=fail );
fi;


#####################################################################
CycleToClass:=function(v)
local u;

if cohdim = 0 then
return [];
fi;

u:=GF2ToZ(SolutionMat(CobandCoc*Z(2),v*Z(2)));
return List([1..cohdim],x->u[Length(u)-cohdim+x]);

end;
#####################################################################

#####################################################################
ClassToCycle:=function(u)
local v,w, i, temp;

w:=List([1..Dimension(n)],x->0);

if cohdim>0 then
for i in [1..Dimension(n)] do
temp := 0;
for j in [1..cohdim] do
temp := temp + CobandCoc[Length(CobandCoc)-cohdim+j][i]*u[j];
w[i] := temp mod 2;
od;
od;
fi;


return w;
end;
#####################################################################

return 	rec(
		#cocyclesBasis:=BasisKerd1,
	 	#boundariesCoefficients:=Rels,
        #torsionCoefficients:=TorsionCoefficients,
        cocyclesBasis:=CobandCoc,
        Mod2Cohomologydim:=cohdim,
	 	cocycleToClass:=CycleToClass,
	 	classToCocycle:=ClassToCycle );

end;
#####################################################################
#####################################################################



#####################################################################
#####################################################################

Mod2CupProduct:=function(arg)
local
    R, u, v, p, q, P, Q, N,
    uCocycle,
    vCocycle,
    uvCocycle,
    uChainMap,
    DimensionR,
    i, w, x, sw;

    ####################BEGIN TO READ THE INPUT##################
R:=arg[1];
DimensionR:=R!.dimension;
u:=arg[2];
v:=arg[3];
p:=arg[4];
q:=arg[5];

if Length(arg)>5 then P:=arg[6];
else
P:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
fi;

if Length(arg)>6 then Q:=arg[7];
else
Q:=CR_Mod2CocyclesAndCoboundaries(R,q,true);
fi;

if Length(arg)>7 then N:=arg[8];
else
N:=CR_Mod2CocyclesAndCoboundaries(R,p+q,true);
fi;
    #####################FINISHED REAQDING THE INPUT#############

uCocycle:=P.classToCocycle(u);
vCocycle:=Q.classToCocycle(v);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);


uvCocycle:=[];
for i in [1..DimensionR(p+q)] do
w:=uChainMap([[i,1]]);
sw:=0;
    for x in w do
    sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
    od;
uvCocycle[i]:=sw mod 2;
od;

return N.cocycleToClass(uvCocycle);
end;
#####################################################################
#####################################################################



#####################################################################
#####################################################################

Mod2RingGenerators:=function(arg)
local
        R, n, GG, IT,
        Gens, GensBasis, Cups, Cupped, cupped, CuppedBasis, spacedim,
        uCocycle, vCocycle, uvCocycle, ww, uChainMap,
        sol, CB, CohomologyBasis, TR,
        BasisP, BasisQ, GF2ToZ,
        i, p, q, u, v, ln, iu, iv, w, x, sw;

n:=arg[2];
spacedim:=3;

if Length(arg)=3 then
spacedim:=arg[3];
fi;

if IsInt(arg[1]) then
IT := arg[1];
if IT = 133 or IT = 138 or IT = 210 or IT = 222 or IT = 224 then
    GG := Image(IsomorphismPcpGroup(SpaceGroupIT(spacedim,IT)));
else
    GG := Image(IsomorphismPcpGroup(SpaceGroupBBNWZ(spacedim,IT)));
fi;
R := ResolutionAlmostCrystalGroup(GG,n+1);

else if IsGroup(arg[1]) then
    GG := Image(IsomorphismPcpGroup(arg[1]));
    R := ResolutionAlmostCrystalGroup(GG,n+1);
    else
        R:=arg[1];
    fi;
fi;


TR:=HomToIntegersModP(R,2);
if Cohomology(TR,n) = 0 then return []; fi;
#####################################################################
CohomologyBasis:=function(Torsion)
local i, v, Basis;
Basis:=[];
for i in [1..Length(Torsion)] do
v:=List([1..Length(Torsion)], j->0);
v[i]:=Torsion[i];
Append(Basis, [v]);
od;
return Basis;
end;
#####################################################################
GF2ToZ:=function(v)
local v0,k;

v0:=[];
for k in [1..Length(v)] do
    if v[k] = 0*Z(2) then
    v0[k]:=0;
    else
        v0[k]:=1;
    fi;
od;

return v0;
end;
#####################################################################

Cups:=CohomologyBasis(List([1..Cohomology(TR,n)],i->1));

if n = 1 then

return Cups;

fi;

CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

Cupped :=[];

for p in [1..QuoInt(n,2)] do
q:=n-p;
BasisP:=CohomologyBasis(List([1..Cohomology(TR,p)],i->1));
BasisQ:=CohomologyBasis(List([1..Cohomology(TR,q)],i->1));


iu :=1;
for u in BasisP do

uCocycle:=CB[p].classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);
ww:=[];
for i in [1..(R!.dimension(n))] do
Append(ww, [uChainMap([[i,1]])]);
od;

iv :=1;
for v in BasisQ do

vCocycle:=CB[q].classToCocycle(v);

if ((p < q) or (p=q and iv>=iu)) then

    uvCocycle:=[];
    for i in [1..(R!.dimension(n))] do
        #w:=uChainMap([[i,1]]);
        w:=ww[i];
        sw:=0;
        for x in w do
            sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;

    cupped := CB[n].cocycleToClass(uvCocycle);
    
    Append(Cupped,[cupped*Z(2)]);

    #cupped :=Mod2CupProduct(R,u,v,p,q,CB[p],CB[q],CB[n]);

fi;

iv := iv+1;
od;

iu := iu+1;
od;

od;

if Cupped = [] then
CuppedBasis := [];
else
CuppedBasis := List(BaseMat(Cupped),ShallowCopy);
fi;

#Append(CuppedBasis,[List([1..Cohomology(TR,n)],x->0*Z(2))]);

#Gens := List(BaseOrthogonalSpaceMat(CuppedBasis),ShallowCopy);

#Print(Cups);

Gens := BaseSteinitzVectors(Cups*Z(2),CuppedBasis)!.factorspace;

GensBasis :=[];

for i in [1..Length(Gens)] do
    GensBasis[i]:=GF2ToZ(Gens[i]);
od;


return GensBasis;
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################

Mod2RingGenerators1to5:=function(arg)
local
        R,n,GG,IT,Gen1,Gen2,Gen3,Gen4,Gen5,Gens,
        cupped, Cups,
        Letters, TR,
        i, p, q, r, u, v, ln, iu, iv, w, x, sw;


if Length(arg)=1 then
n:=5;
else
n:=arg[2];
fi;


if IsInt(arg[1]) then
IT := arg[1];
if IT = 133 or IT = 138 or IT = 210 or IT = 222 or IT = 224 then
    GG := Image(IsomorphismPcpGroup(SpaceGroupIT(3,IT)));
else
    GG := Image(IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT)));
fi;
R := ResolutionAlmostCrystalGroup(GG,n+1);

else if IsGroup(arg[1]) then
    GG := Image(IsomorphismPcpGroup(arg[1]));
    R := ResolutionAlmostCrystalGroup(GG,n+1);
    else
        R:=arg[1];
    fi;
fi;

TR:=HomToIntegersModP(R,2);

#####################################################################



Letters := [["A1","A2","A3","A4","A5","A6","A7"],["B1","B2","B3","B4","B5","B6","B7"],["C1","C2","C3","C4","C5","C6","C7"],["D1","D2","D3","D4","D5","D6","D7"]];

Gen1:=Mod2RingGenerators(R,1);
Gen2:=Mod2RingGenerators(R,2);
Gen3:=Mod2RingGenerators(R,3);
Gen4:=Mod2RingGenerators(R,4);
Gen5:=Mod2RingGenerators(R,5);

Gens:=[Gen1,Gen2,Gen3,Gen4,Gen5];


return [IT, Gens];
end;

#####################################################################
#####################################################################





#####################################################################
#####################################################################


Mod2RingGensAndRels:=function(arg)
local
        R,n,GG,IT,Gen1,Gen2,Gen3,Gen4,spacedim,
        Gens, Cupped, CuppedLetter, CupRels, CupRelsLetter,
        CupBase2, CupBase3,CupBase4,CupBase5,CupBase6,
        CupBase2Letter,CupBase3Letter,CupBase4Letter,CupBase5Letter,CupBase6Letter,
        CupRel2, CupRel2Letter, CupRel3Letter, CupRel4Letter, CupRel5Letter, CupRel6Letter,
        cupped, Cups,
        Letters, AddLetters,
        uCocycle, vCocycle, uvCocycle, uChainMap, ww,
        sol, cc, CB, CohomologyBasis, TR,
        BasisP, BasisQ, SmithRecord, GF2ToZ, IToPosition,
        i, p, q, r, u, v, ln, iu, iv, w, x, sw;


if Length(arg)=1 then
n:=6;
spacedim:=3;
else
n:=arg[2];
spacedim:=3;
fi;

if Length(arg)=3 then
n:=6;
spacedim:=arg[3];
fi;

if IsInt(arg[1]) then
IT := arg[1];
if IT = 133 or IT = 138 or IT = 210 or IT = 222 or IT = 224 then
    GG := Image(IsomorphismPcpGroup(SpaceGroupIT(spacedim,IT)));
else
    GG := Image(IsomorphismPcpGroup(SpaceGroupBBNWZ(spacedim,IT)));
fi;
R := ResolutionAlmostCrystalGroup(GG,n+1);

else if IsGroup(arg[1]) then
    GG := Image(IsomorphismPcpGroup(arg[1]));
    R := ResolutionAlmostCrystalGroup(GG,n+1);
    else
        R:=arg[1];
    fi;
fi;

TR:=HomToIntegersModP(R,2);

#####################################################################
CohomologyBasis:=function(Torsion)
local i, v, Basis;
Basis:=[];
for i in [1..Length(Torsion)] do
v:=List([1..Length(Torsion)], j->0);
v[i]:=Torsion[i];
Append(Basis, [v]);
od;
return Basis;
end;
#####################################################################
GF2ToZ:=function(v)
local v0,k;

v0:=[];
for k in [1..Length(v)] do
    if v[k] = 0*Z(2) then
    v0[k]:=0;
    else
        v0[k]:=1;
    fi;
od;

return v0;
end;
######################################################################
IToPosition:=function(v)
local v0,k;

v0:=[];
for k in [1..Length(v)] do
    if v[k] = 1 then
    Append(v0, [ k ]);
    fi;
od;
return v0;
end;
######################################################################
AddLetters:=function(v)
local v0,k;

if v = [] then
 return [];
fi;
v0:=v[1];
if Length(v)>1 then
for k in [2..Length(v)] do
v0:=Concatenation(v0,"+",v[k]);
od;
fi;
return [v0];
end;
######################################################################

Letters := [["A1","A2","A3","A4","A5","A6","A7"],["B1","B2","B3","B4","B5","B6","B7"],["C1","C2","C3","C4","C5","C6","C7"],["D1","D2","D3","D4","D5","D6","D7"]];

Gen1:=Mod2RingGenerators(R,1,spacedim);
Gen2:=Mod2RingGenerators(R,2,spacedim);
Gen3:=Mod2RingGenerators(R,3,spacedim);
Gen4:=Mod2RingGenerators(R,4,spacedim);
Gens:=[Gen1,Gen2,Gen3,Gen4];

CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;


####################### r = 2 ##########################

CupBase2 := [];
CupBase2Letter := [];

iu :=1;
for u in Gen1 do
iv :=1;
for v in Gen1 do
if iv>=iu then
    cupped :=Mod2CupProduct(R,u,v,1,1,CB[1],CB[1],CB[2]);
    Append(CupBase2,[cupped]);
    Append(CupBase2Letter,[Concatenation(Letters[1,iu],Letters[1,iv])]);
fi;
iv := iv+1;
od;
iu := iu+1;
od;


CupRel2 := [];
CupRel2Letter := [];

if not (CupBase2 = []) then
CupRel2 := List(BasisNullspaceModN(CupBase2,2),ShallowCopy)*Z(2);
for cc in CupRel2 do
    Append(CupRel2Letter,AddLetters(List(IToPosition(GF2ToZ(cc)),x->CupBase2Letter[x])));
od;
fi;


#CupBase2 := List(BaseOrthogonalSpaceMat(CupRels),ShallowCopy);
#CupBase2 := List(BaseMat(TransposedMat(Cupped)*Z(2)),ShallowCopy);
####Both are problematic!!!

#iu :=1;
#for cc in Gen2 do
#    Append(CupBase2, [cc]);
#    Append(CupBase2Letter,[Letters[2,iu]]);
#    iu :=iu+1;
#od;



####################### r = 3 ##########################



CupBase3 :=[];
CupBase3Letter :=[];
iu :=1;
for u in Gen1 do
iv :=1;
for v in CupBase2 do
    cupped :=Mod2CupProduct(R,u,v,1,2,CB[1],CB[2],CB[3]);
    Append(CupBase3,[cupped]);
    Append(CupBase3Letter,[Concatenation(Letters[1,iu],CupBase2Letter[iv])]);
    iv := iv+1;
od;
iu := iu+1;
od;


CupRel3Letter := [];

iu :=1;
for u in Gen1 do
iv :=1;
for cc in Gen2 do
    cupped :=Mod2CupProduct(R,u,cc,1,2,CB[1],CB[2],CB[3]);
    if cupped = List([1..Cohomology(TR,3)],x->0) then
        Append(CupRel3Letter,[Concatenation(Letters[1,iu],Letters[2,iv])]);
    else
    if not (CupBase3 = []) then
    sol :=SolutionMat(CupBase3*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase3,[cupped]);
        Append(CupBase3Letter,[Concatenation(Letters[1,iu],Letters[2,iv])]);
    else
        #Print(sol);
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase3Letter[x]);
        Append(sol,[Concatenation(Letters[1,iu],Letters[2,iv])]);
        Append(CupRel3Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;



####################### r = 4 ##########################

CupBase4 :=[];
CupBase4Letter :=[];
iu :=1;
for u in Gen1 do
iv :=1;
for v in CupBase3 do
    cupped :=Mod2CupProduct(R,u,v,1,3,CB[1],CB[3],CB[4]);
    Append(CupBase4,[cupped]);
    Append(CupBase4Letter,[Concatenation(Letters[1,iu],CupBase3Letter[iv])]);
    iv := iv+1;
od;
iu := iu+1;
od;


CupRel4Letter := [];



iu :=1;
for u in Gen1 do
iv :=1;
for cc in Gen3 do
    cupped :=Mod2CupProduct(R,u,cc,1,3,CB[1],CB[3],CB[4]);
    if cupped = List([1..Cohomology(TR,4)],x->0) then
        Append(CupRel4Letter,[Concatenation(Letters[1,iu],Letters[3,iv])]);
    else
    if not (CupBase4 = []) then
    sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase4,[cupped]);
        Append(CupBase4Letter,[Concatenation(Letters[1,iu],Letters[3,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase4Letter[x]);
        Append(sol,[Concatenation(Letters[1,iu],Letters[3,iv])]);
        Append(CupRel4Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;


    
iu :=1;
for u in Gen2 do
iv :=1;
for cc in Gen2 do
    if iv>=iu then
    cupped :=Mod2CupProduct(R,u,cc,2,2,CB[2],CB[2],CB[4]);
    if cupped = List([1..Cohomology(TR,4)],x->0) then
        Append(CupRel4Letter,[Concatenation(Letters[2,iu],Letters[2,iv])]);
    else
    if not (CupBase4 = []) then
    sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase4,[cupped]);
        Append(CupBase4Letter,[Concatenation(Letters[2,iu],Letters[2,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase4Letter[x]);
        Append(sol,[Concatenation(Letters[2,iu],Letters[2,iv])]);
        Append(CupRel4Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;


####################### r = 5 ##########################

CupBase5 :=[];
CupBase5Letter :=[];
iu :=1;
for u in Gen1 do
iv :=1;
for v in CupBase4 do
    cupped :=Mod2CupProduct(R,u,v,1,4,CB[1],CB[4],CB[5]);
    Append(CupBase5,[cupped]);
    Append(CupBase5Letter,[Concatenation(Letters[1,iu],CupBase4Letter[iv])]);
    iv := iv+1;
od;
iu := iu+1;
od;


CupRel5Letter := [];



iu :=1;
for u in Gen1 do

uCocycle:=CB[1].classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,1,4);
ww:=[];
for i in [1..(R!.dimension(5))] do
Append(ww, [uChainMap([[i,1]])]);
od;


iv :=1;
for cc in Gen4 do
    
    vCocycle:=CB[4].classToCocycle(cc);

    uvCocycle:=[];
    for i in [1..(R!.dimension(5))] do
        w:=ww[i];
        sw:=0;
        for x in w do
            sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;

    cupped := CB[5].cocycleToClass(uvCocycle);
    #cupped :=Mod2CupProduct(R,u,cc,1,4,CB[1],CB[4],CB[5]);
    if cupped = List([1..Cohomology(TR,5)],x->0) then
        Append(CupRel5Letter,[Concatenation(Letters[1,iu],Letters[4,iv])]);
    else
    if not (CupBase5 = []) then
    sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase5,[cupped]);
        Append(CupBase5Letter,[Concatenation(Letters[1,iu],Letters[4,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase5Letter[x]);
        Append(sol,[Concatenation(Letters[1,iu],Letters[4,iv])]);
        Append(CupRel5Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;

iu :=1;
for u in Gen2 do

uCocycle:=CB[2].classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,2,3);
ww:=[];
for i in [1..(R!.dimension(5))] do
Append(ww, [uChainMap([[i,1]])]);
od;


iv :=1;
for cc in Gen3 do

    vCocycle:=CB[3].classToCocycle(cc);

    uvCocycle:=[];
    for i in [1..(R!.dimension(5))] do
        w:=ww[i];
        sw:=0;
        for x in w do
            sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;

    cupped := CB[5].cocycleToClass(uvCocycle);
    #cupped :=Mod2CupProduct(R,u,cc,2,3,CB[2],CB[3],CB[5]);

    if cupped = List([1..Cohomology(TR,5)],x->0) then
        Append(CupRel5Letter,[Concatenation(Letters[2,iu],Letters[3,iv])]);
    else
    if not (CupBase5 = []) then
    sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase5,[cupped]);
        Append(CupBase5Letter,[Concatenation(Letters[2,iu],Letters[3,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase5Letter[x]);
        Append(sol,[Concatenation(Letters[2,iu],Letters[3,iv])]);
        Append(CupRel5Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;


####################### r = 6 ##########################

CupBase6 :=[];
CupBase6Letter :=[];
iu :=1;
for u in Gen1 do
iv :=1;
for v in CupBase5 do
    cupped :=Mod2CupProduct(R,u,v,1,5,CB[1],CB[5],CB[6]);
    Append(CupBase6,[cupped]);
    Append(CupBase6Letter,[Concatenation(Letters[1,iu],CupBase5Letter[iv])]);
    iv := iv+1;
od;
iu := iu+1;
od;

iu :=1;
for u in Gen2 do
iv :=1;
for v in CupBase4 do
    cupped :=Mod2CupProduct(R,u,v,2,4,CB[2],CB[4],CB[6]);
    Append(CupBase6,[cupped]);
    Append(CupBase6Letter,[Concatenation(Letters[2,iu],CupBase4Letter[iv])]);
    iv := iv+1;
od;
iu := iu+1;
od;


CupRel6Letter := [];



iu :=1;
for u in Gen2 do

uCocycle:=CB[2].classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,2,4);
ww:=[];
for i in [1..(R!.dimension(6))] do
Append(ww, [uChainMap([[i,1]])]);
od;


iv :=1;
for cc in Gen4 do

    vCocycle:=CB[4].classToCocycle(cc);

    uvCocycle:=[];
    for i in [1..(R!.dimension(6))] do
        w:=ww[i];
        sw:=0;
        for x in w do
            sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;

    cupped := CB[6].cocycleToClass(uvCocycle);
    #cupped :=Mod2CupProduct(R,u,cc,2,4,CB[2],CB[4],CB[6]);

    if cupped = List([1..Cohomology(TR,6)],x->0) then
        Append(CupRel6Letter,[Concatenation(Letters[2,iu],Letters[4,iv])]);
    else
    if not (CupBase6 = []) then
    sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase6,[cupped]);
        Append(CupBase6Letter,[Concatenation(Letters[2,iu],Letters[4,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase6Letter[x]);
        Append(sol,[Concatenation(Letters[2,iu],Letters[4,iv])]);
        Append(CupRel6Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;

iu :=1;
for u in Gen3 do

uCocycle:=CB[3].classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,3);
ww:=[];
for i in [1..(R!.dimension(6))] do
Append(ww, [uChainMap([[i,1]])]);
od;

iv :=1;
for cc in Gen3 do
    if iv>=iu then
    vCocycle:=CB[3].classToCocycle(cc);

    uvCocycle:=[];
    for i in [1..(R!.dimension(6))] do
        w:=ww[i];
        sw:=0;
        for x in w do
            sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;

    cupped := CB[6].cocycleToClass(uvCocycle);
    #cupped :=Mod2CupProduct(R,u,cc,3,3,CB[3],CB[3],CB[6]);
        
    if cupped = List([1..Cohomology(TR,6)],x->0) then
        Append(CupRel6Letter,[Concatenation(Letters[3,iu],Letters[3,iv])]);
    else
    if not (CupBase6 = []) then
    sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
    if sol = fail then
        Append(CupBase6,[cupped]);
        Append(CupBase6Letter,[Concatenation(Letters[3,iu],Letters[3,iv])]);
    else
        sol:=List(IToPosition(GF2ToZ(sol)),x->CupBase6Letter[x]);
        Append(sol,[Concatenation(Letters[3,iu],Letters[3,iv])]);
        Append(CupRel6Letter,AddLetters(sol));
    fi;
    fi;
    fi;
    fi;
    iv := iv+1;
od;
iu := iu+1;
od;




return rec(GensAtDegN:=List([1..4],x->List([1..Length(Gens[x])],y->Letters[x,y])),RelsAtDegN:=[CupRel2Letter,CupRel3Letter,CupRel4Letter,CupRel5Letter,CupRel6Letter]);
end;
#####################################################################
#####################################################################
