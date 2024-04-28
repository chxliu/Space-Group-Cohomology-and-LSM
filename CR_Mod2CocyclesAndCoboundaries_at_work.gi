Read("~/Downloads/Space_Group_Cocycles.gi");



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
    #sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
        sw:=sw+ vCocycle[AbsoluteValue(x[1])];
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
        Gens, GensBasis, GensBasis1ton, Cups, Cupped, cupped, CuppedBasis, spacedim,
        uCocycle, vCocycle, uvCocycle, ww, uChainMap,
        sol, CB, CohomologyBasis, TR,
        BasisP, BasisQ, GF2ToZ,
        i, j, p, q, u, v, ln, iu, iv, w, x, sw;

#This function computes, for a given n, the generators at degree 1,2,...,n.
#e.g. Mod2RingGenerators(IT=76,n=4);

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

else
    if IsGroup(arg[1]) then
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
#####################################################################


CB:=[];
for p in [1..n] do
    CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);  #CR_Mod2CocyclesAndCoboundaries gives all the cocycles followed by all the coboundaries as vectors
od;


GensBasis1ton :=[CohomologyBasis(List([1..Cohomology(TR,1)],i->1))];  # Record degree-1 generators  #Make an identity matrix, consisting of basis vectors for cocycles


for j in [2..n] do        # Then deal with degree-j generators for j=2,3,...,n

    Cups:=CohomologyBasis(List([1..Cohomology(TR,j)],i->1));    #Make an identity matrix, consisting of basis vectors for cocycles

    Cupped :=[];

    for p in [QuoInt(j+1,2)..(j-1)] do  #QuoInt(j,2) is gap command for the usual Int(j/2)
        q:=j-p;
        BasisP:=CohomologyBasis(List([1..Cohomology(TR,p)],i->1));
        BasisQ:=CohomologyBasis(List([1..Cohomology(TR,q)],i->1));

        iu :=1;
        for u in BasisP do

            uCocycle:=CB[p].classToCocycle(u);
            uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);
            ww:=[];
            for i in [1..(R!.dimension(j))] do
                Append(ww, [uChainMap([[i,1]])]);
            od;

            iv :=1;
            for v in BasisQ do

                vCocycle:=CB[q].classToCocycle(v);

                if ((p > q) or (p=q and iv>=iu)) then

                    uvCocycle:=[];
                    for i in [1..(R!.dimension(j))] do
                        w:=ww[i];
                        sw:=0;
                        for x in w do
                            sw:=sw + vCocycle[AbsoluteValue(x[1])];
                        od;
                        uvCocycle[i]:=sw mod 2;
                    od;

                    cupped := CB[j].cocycleToClass(uvCocycle);
    
                    Append(Cupped,[cupped*Z(2)]);

                    #cupped :=Mod2CupProduct(R,u,v,p,q,CB[p],CB[q],CB[j]);

                fi;

                iv := iv+1;
            od;

            iu := iu+1;
        od;

    od;

    if Cups = [] then        #This is when cohomology dim = 0
        GensBasis1ton[j]:= [];
    else
        CuppedBasis := List(BaseMat(Cupped),ShallowCopy);
        Gens := BaseSteinitzVectors(Cups*Z(2),CuppedBasis)!.factorspace;

        GensBasis :=[];

        for i in [1..Length(Gens)] do
            GensBasis[i]:=GF2ToZ(Gens[i]);
        od;

        GensBasis1ton[j]:=GensBasis;
    fi;

od;

return GensBasis1ton;
end;

#####################################################################
#####################################################################



#####################################################################
#####################################################################


Mod2RingGensAndRels:=function(arg)
local
        R,n,GG,IT,Gen1,Gen2,Gen3,Gen4,Gen5,spacedim,GenDim1to4,
        Gens, GensLett, Cupped, CupRelsLett, CupTemp, CupTempLett,
        CupBase2all, CupBase2, CupBase3,CupBase4,CupBase5,CupBase6,
        CupBase2Lett,CupBase3Lett,CupBase4Lett,CupBase5Lett,CupBase6Lett,
        CupRel2Lett, CupRel3Lett, CupRel4Lett, CupRel5Lett, CupRel6Lett,
        cupped,
        Lett1, Lett2, IO,
        uCocycle, vCocycle, uvCocycle, uChainMap, ww,
        sol, solrel, cc, CB, CohomologyBasis, TR,
        BasisP, BasisQ, SmithRecord, GF2ToZ, IToPosition,
        NonNegativeVec, Letter2Monomial, PrintMonomialString,
        i, p, q, r, u, v, ln, iu, iv, w, x, sw;

#Standard input: arg[1] = IT (# of space group), arg[2] = n (relations up to deg(n) is calculated)
#e.g.: Mod2RingGensAndRels(89);
#e.g.: Mod2RingGensAndRels(89,3);
#e.g.: Mod2RingGensAndRels(89,3,R89);
#e.g.: Mod2RingGensAndRels(89,3,R89,Gens);


n:=6; #This is the highest degree at which the relations are calculated. For now we fix it to be 6.

if Length(arg)=1 then
spacedim:=3;
fi;

if Length(arg)>=2 then
spacedim:=arg[2];
fi;


IT := arg[1];  #Group number in International Table for Crystallography

Print("===========================================\n");
Print("Begin Group No. ", IT, ":\n");

if IT = 219 then
    Print("Caution: the group being calculated (No. 219) has two degree-6 generators!\n");
    Print("These degree-6 generators have not been included in the current code, ");
    Print("Current program only outputs relations at degree <=6. Higher (>6) degree relations must be calculated manually.\n");
fi;
if IT = 226 or IT = 228 then
    Print("Caution: the group being calculated (either No. 226 or No. 228) has a degree-6 generator!\n");
    Print("This degree-6 generator has not been included in the current code, ");
    Print("as the length-7 resolution required to get the degree-6 generator exceeds the current memory limit.\n");
    Print("Current program only outputs relations at degree <=6. Higher (>6) degree relations must be calculated manually.\n");
fi;




if Length(arg)<=2 then

    if IT = 133 or IT = 138 or IT = 210 or IT = 222 or IT = 224 then
        GG := Image(IsomorphismPcpGroup(SpaceGroupIT(spacedim,IT)));
    else
        GG := Image(IsomorphismPcpGroup(SpaceGroupBBNWZ(spacedim,IT)));
    fi;
    
    R := ResolutionAlmostCrystalGroup(GG,n+1);            #Construct resolution;
    Gens:=Mod2RingGenerators(R,4,spacedim);               #Calculate generators
    

else                                                      #otherwise, Length(arg) = 4

    R := arg[3];                                          #receive resolution from input;
    Gens:=arg[4];                                         #receive generators from input;

fi;


Gen1:=Gens[1];
Gen2:=Gens[2];
Gen3:=Gens[3];
Gen4:=Gens[4];
Gen5:=[]; #based on the a posteriori fact that no space group has degree-5 generators.
GenDim1to4:=[Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4)];

if Length(Gen4)>0 then
    Print("Caution: this group contains degree-4 generators!\n");
    Print("Current program only outputs relations at degree <=6. Higher (>6) degree relations must be calculated manually.\n");
fi;


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
NonNegativeVec:=function(v)
local ps,k;

ps := true;
for k in [1..Length(v)] do
    if v[k] < 0 then
    ps := false;
    fi;
od;
return ps;
end;
######################################################################
Letter2Monomial:=function(vec,GensDim)
local j,k,u,v,lett;

lett := [];
u := 1;
for j in [1..Length(GensDim)] do
    for k in [1..GensDim[j]] do
        if vec[u]>1 then
            Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k),"^",vec[u]],"")]);
        else
            if vec[u]=1 then
                Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k)],"")]);
            fi;
        fi;
        u := u+1;
    od;
od;
return JoinStringsWithSeparator(lett, ".");
end;
######################################################################
PrintMonomialString:=function(vecs,GensDim,sep)
local poly;

poly:=List(vecs,x->Letter2Monomial(x,GensDim));
Print(JoinStringsWithSeparator(poly,sep)," ");
return 0;
end;
######################################################################




GensLett:=CohomologyBasis(List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))],i->1));
#GensLett records the powers, later to be used to convert generators/relations to letters.


TR:=HomToIntegersModP(R,2);            #Apply the Hom functor.


CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

Print("Number of Generators at degrees 1-4: ");
Print(GenDim1to4);
Print("\n");


####################### r = 1 ##########################

#Print("Chosen basis at degree 1:\n");
#PrintMonomialString(List([1..Length(Gen1)],x->GensLett[x]),GenDim1to4,",");

Print("Matching dimensions... [dim(Chosen basis) = dim(H)]\n");
Print("dim(H^1)=", Cohomology(TR,1),", ");

####################### r = 2 ##########################

CupRelsLett := [];

CupBase2 := [];
CupBase2Lett := [];
CupRel2Lett := [];

iu :=1;                                           #First step: calculate u-cup-v for u != v
for u in Gen1 do
    iv :=1;
    for v in Gen1 do
        Lett1 := GensLett[iu]+GensLett[iv];
        if iv>iu then
            cupped := Mod2CupProduct(R,u,v,1,1,CB[1],CB[1],CB[2]);
    
            if cupped = List([1..Cohomology(TR,2)],x->0) then   #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel2Lett,[[Lett1]]);
            else
                if CupBase2 = [] then             #if no basis yet then push in the genuine cocycle u-cup-v
                    Append(CupBase2,[cupped]);
                    Append(CupBase2Lett,[Lett1]);
                else
                    sol :=SolutionMat(CupBase2*Z(2),cupped*Z(2));
                    if sol = fail then                        #if u-cup-v is a genuine new cocycle
                        Append(CupBase2,[cupped]);
                        Append(CupBase2Lett,[Lett1]);
                    else                       #if u-cup-v is expressable by other cocycles
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                        if (Lett1 in CupBase2Lett) = false then
                            Append(CupRelsLett,[Lett1]);
                            Append(CupRel2Lett,[solrel]);
                        fi;
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;


iu :=1;                                    #Second step: calculate u-cup-u
for u in Gen1 do
    Lett1 := GensLett[iu]+GensLett[iu];
    cupped := Mod2CupProduct(R,u,u,1,1,CB[1],CB[1],CB[2]);
    
        if cupped = List([1..Cohomology(TR,2)],x->0) then   #if u-cup-u is a coboundary
            Append(CupRelsLett,[Lett1]);
            Append(CupRel2Lett,[[Lett1]]);
        else
            if CupBase2 = [] then             #if no basis yet then push in the genuine cocycle u-cup-u
                Append(CupBase2,[cupped]);
                Append(CupBase2Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase2*Z(2),cupped*Z(2));
                if sol = fail then                                  #if u-cup-u is a genuine new cocycle
                    Append(CupBase2,[cupped]);
                    Append(CupBase2Lett,[Lett1]);
                else                              #if u-cup-u is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                    if (Lett1 in CupBase2Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel2Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
    iu := iu+1;
od;

Append(CupBase2,Gen2);
for iu in [(1+Length(Gen1))..(Length(Gen1)+Length(Gen2))] do
    Append(CupBase2Lett,[GensLett[iu]]);
od;


#Print("Independent relations at degree 2:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 2:\n");
#PrintMonomialString(CupBase2Lett,GenDim1to4,",");

if Length(CupBase2Lett) = Cohomology(TR,2) then
    Print("dim(H^2)=", Cohomology(TR,2),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^2) = ", Length(CupBase2Lett) - Cohomology(TR,2),"\n");
fi;




####################### r = 3 ##########################


CupBase3 :=[];
CupBase3Lett :=[];
CupRel3Lett := [];
CupTemp := [];
CupTempLett := [];

#Begins: degree-2 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase2 do
    iv :=1;
    for v in Gen1 do
        
        cupped :=Mod2CupProduct(R,u,v,2,1,CB[2],CB[1],CB[3]);      #then calculate u-cup-v
        
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        Lett1 := CupBase2Lett[iu] + GensLett[iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;                      #then rewrite IO
                Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                Append(CupTempLett,[Lett1]);
            fi;
        od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations

        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
            if cupped = List([1..Cohomology(TR,3)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel3Lett,[[Lett1]]);
            elif CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase3,[cupped]);
                Append(CupBase3Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase3*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase3,[cupped]);
                    Append(CupBase3Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase3Lett[x]));
                    if (Lett1 in CupBase3Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel3Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Finished: degree-2 cup with degree-1-gen


Append(CupBase3,Gen3);
for iu in [(1+Length(Gen1)+Length(Gen2))..(Length(Gen1)+Length(Gen2)+Length(Gen3))] do
    Append(CupBase3Lett,[GensLett[iu]]);
od;

if CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
    if (CupTemp = []) = false then
        Append(CupBase3,[CupTemp[1]]);
        Append(CupBase3Lett,[CupTempLett[1]]);
    fi;
fi;

if (CupBase3 = []) = false then

    for cupped in CupTemp do

        sol :=SolutionMat(CupBase3*Z(2),cupped*Z(2));
        if sol = fail then                  #if u-cup-v is a genuine new cocycle
            Append(CupBase3,[cupped]);
            Append(CupBase3Lett,[Lett1]);
        fi;
    od;
fi;

#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 3:\n");
#PrintMonomialString(CupBase3Lett,GenDim1to4,",");

if Length(CupBase3Lett) = Cohomology(TR,3) then
    Print("dim(H^3)=", Cohomology(TR,3),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^3) = ", Length(CupBase3Lett) - Cohomology(TR,3),"\n");
fi;




####################### r = 4 ##########################

CupBase4 :=[];
CupBase4Lett :=[];
CupRel4Lett := [];
CupTemp := [];
CupTempLett := [];

#Step-1 begins here: degree-3 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase3 do
    iv :=1;
    for v in Gen1 do

        cupped :=Mod2CupProduct(R,u,v,3,1,CB[3],CB[1],CB[4]);      #then calculate u-cup-v
    
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        Lett1 := CupBase3Lett[iu] + GensLett[iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;                      #then rewrite IO
                Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                Append(CupTempLett,[Lett1]);
            fi;
        od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations
    
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,

            if cupped = List([1..Cohomology(TR,4)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel4Lett,[[Lett1]]);
            elif CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                    Append(CupBase4,[cupped]);
                    Append(CupBase4Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase4,[cupped]);
                    Append(CupBase4Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
                    if (Lett1 in CupBase4Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel4Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-3 cup with degree-1-gen



if CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
    if (CupTemp = []) = false then
        Append(CupBase4,[CupTemp[1]]);
        Append(CupBase4Lett,[CupTempLett[1]]);
    fi;
fi;

if (CupBase4 = []) = false then

    for cupped in CupTemp do

        sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
        if sol = fail then                  #if u-cup-v is a genuine new cocycle
            Append(CupBase4,[cupped]);
            Append(CupBase4Lett,[Lett1]);
        fi;
    od;
fi;



#Step-2 begins here: degree-2-gen cup with degree-2-gen
#
#
iu :=1;
for u in Gen2 do
    iv :=1;
    for v in Gen2 do
        if iv>=iu then
        
            cupped :=Mod2CupProduct(R,u,v,2,2,CB[2],CB[2],CB[4]);      #then calculate u-cup-v
        
            Lett1 := GensLett[Length(Gen1)+iu] + GensLett[Length(Gen1)+iv];
        
            if cupped = List([1..Cohomology(TR,4)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel4Lett,[[Lett1]]);
            elif CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase4,[cupped]);
                Append(CupBase4Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase4,[cupped]);
                    Append(CupBase4Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
                    if (Lett1 in CupBase4Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel4Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-2 finished: degree-2-gen cup with degree-2-gen

Append(CupBase4,Gen4);
for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))] do
    Append(CupBase4Lett,[GensLett[iu]]);
od;




#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 4:\n");
#PrintMonomialString(CupBase4Lett,GenDim1to4,",");

if Length(CupBase4Lett) = Cohomology(TR,4) then
    Print("dim(H^4)=", Cohomology(TR,4),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^4) = ", Length(CupBase4Lett) - Cohomology(TR,4),"\n");
fi;



####################### r = 5 ##########################

CupBase5 :=[];
CupBase5Lett :=[];
CupRel5Lett := [];
CupTemp := [];
CupTempLett := [];

#Step-1 begins here: degree-4 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase4 do

    #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,1);
    ww:=[];
    for i in [1..(R!.dimension(5))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        #cupped :=Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]);      #then calculate u-cup-v
        
        #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(5))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[5].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 2 ended ####


        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        Lett1 := CupBase4Lett[iu] + GensLett[iv];
        IO := false;
        
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;                      #then rewrite IO
                Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                Append(CupTempLett,[Lett1]);
            fi;
        od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations

        
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        
            if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel5Lett,[[Lett1]]);
            elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase5,[cupped]);
                    Append(CupBase5Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
                    if (Lett1 in CupBase5Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel5Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-4 cup with degree-1-gen


#Step-2 begins here: degree-3-gen cup with degree-2-gen
#
#
iu :=1;
for u in Gen3 do

    #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 1: ####
    ####
    uCocycle:=CB[3].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,2);
    ww:=[];
    for i in [1..(R!.dimension(5))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5])  -- part 1 ended ####

    iv :=1;
    for v in Gen2 do
        
        #cupped :=Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]);      #then calculate u-cup-v
            
        #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 2: ####
        ####
        vCocycle:=CB[2].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(5))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[5].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 2 ended ####
        
        
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;                      #then rewrite IO
                Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                Append(CupTempLett,[Lett1]);
            fi;
        od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations

        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,

            if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel5Lett,[[Lett1]]);
            elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase5,[cupped]);
                    Append(CupBase5Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
                    if (Lett1 in CupBase5Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel5Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-2 finished: degree-3-gen cup with degree-2-gen


Append(CupBase5,Gen5);
for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))] do
    Append(CupBase5Lett,[GensLett[iu]]);
od;


if CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
    if (CupTemp = []) = false then
        Append(CupBase5,[CupTemp[1]]);
        Append(CupBase5Lett,[CupTempLett[1]]);
    fi;
fi;

if (CupBase5 = []) = false then

    for cupped in CupTemp do

        sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
        if sol = fail then                  #if u-cup-v is a genuine new cocycle
            Append(CupBase5,[cupped]);
            Append(CupBase5Lett,[Lett1]);
        fi;
    od;
fi;


#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 5:\n");
#PrintMonomialString(CupBase5Lett,GenDim1to4,",");


if Length(CupBase5Lett) = Cohomology(TR,5) then
        Print("dim(H^5)=", Cohomology(TR,5),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^5) = ", Length(CupBase5Lett) - Cohomology(TR,5),"\n");
fi;



####################### r = 6 ##########################

CupBase6 :=[];
CupBase6Lett :=[];
CupRel6Lett := [];
CupTemp := [];
CupTempLett := [];

#Step-1 begins here: degree-5 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase5 do

    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[5].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,5,1);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        #cupped :=Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]);      #then calculate u-cup-v
        
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(6))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[6].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2 ended ####
        
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        Lett1 := CupBase5Lett[iu] + GensLett[iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;                      #then rewrite IO
                Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                Append(CupTempLett,[Lett1]);
            fi;
        od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations
        
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
                
            if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                    if (Lett1 in CupBase6Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel6Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-5 cup with degree-1-gen


#Step-2 begins here: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen
#
#
iu :=1;
for u in CupBase4 do
    if List([1..Length(Gen1)],x->CupBase4Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Bi-cup-Bj or D
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 1: ####
        ####
        uCocycle:=CB[4].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,2);
        ww:=[];
        for i in [1..(R!.dimension(6))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
            #cupped :=Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]);      #then calculate u-cup-v
        
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[6].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2 ended ####
        
            
            #### begin checking if u-cup-v contains letters of the lower-degree relations
            ####
            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];
            IO := false;
            for Lett2 in CupRelsLett do
                if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                    IO :=true;                      #then rewrite IO
                    Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
                    Append(CupTempLett,[Lett1]);
                fi;
            od;
            ####
            #### finished checking if u-cup-v contains letters of the lower-degree relations
            
            
            if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        
                if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel6Lett,[[Lett1]]);
                elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else
                    sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBase6,[cupped]);
                        Append(CupBase6Lett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                        if (Lett1 in CupBase6Lett) = false then
                            Append(CupRelsLett,[Lett1]);
                            Append(CupRel6Lett,[solrel]);
                        fi;
                    fi;
                fi;
            fi;
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen




if CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
    if (CupTemp = []) = false then
        Append(CupBase6,[CupTemp[1]]);
        Append(CupBase6Lett,[CupTempLett[1]]);
    fi;
fi;

if (CupBase6 = []) = false then

    for cupped in CupTemp do

        sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
        if sol = fail then                  #if u-cup-v is a genuine new cocycle
            Append(CupBase6,[cupped]);
            Append(CupBase6Lett,[Lett1]);
        fi;
    od;
fi;



#Step-3 begins here: degree-3-gen cup with degree-3-gen
#
#
iu :=1;
for u in Gen3 do

    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[3].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,3);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen3 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+Length(Gen2)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
            #IO := false;
            #for Lett2 in CupRelsLett do
            #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            #    IO :=true;       #then rewrite IO
            #    fi;
            #od;
        
            #cupped :=Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]);      #then calculate u-cup-v
            
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[3].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[6].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2 ended ####
        
            if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                    if (Lett1 in CupBase6Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel6Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-3-gen cup with degree-3-gen

#Append(CupBase6,Gen6);
#for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))] do
#    Append(CupBase6Lett,[GensLett[iu]]);
#od;



#Print("Chosen basis at degree 6:\n");
#PrintMonomialString(CupBase6Lett,GenDim1to4,",");

if Length(CupBase6Lett) = Cohomology(TR,6) then
        Print("dim(H^6)=", Cohomology(TR,6),".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^6) = ", Length(CupBase6Lett) - Cohomology(TR,6),"\n");
fi;


Print("Independent relations:\n");
#Print(CupRelsLett,"\n");
Print(Length(CupRel2Lett)," at deg 2: ");List(CupRel2Lett,x->PrintMonomialString(x,GenDim1to4,"+"));Print("\n");
Print(Length(CupRel3Lett)," at deg 3: ");List(CupRel3Lett,x->PrintMonomialString(x,GenDim1to4,"+"));Print("\n");
Print(Length(CupRel4Lett)," at deg 4: ");List(CupRel4Lett,x->PrintMonomialString(x,GenDim1to4,"+"));Print("\n");
Print(Length(CupRel5Lett)," at deg 5: ");List(CupRel5Lett,x->PrintMonomialString(x,GenDim1to4,"+"));Print("\n");
Print(Length(CupRel6Lett)," at deg 6: ");List(CupRel6Lett,x->PrintMonomialString(x,GenDim1to4,"+"));Print("\n");


Print("End Group No. ", IT, ".\n");
Print("===========================================\n");

return true;
#rec(GensAtDegN:=List([1..4],x->List([1..Length(Gens[x])],y->Letters[x,y])),RelsAtDegN:=[CupRel2Lett,CupRel3Lett,CupRel4Lett,CupRel5Lett,CupRel6Lett]);
end;
#####################################################################
#####################################################################







#####################################################################
#####################################################################

SpaceGroupCohomologyRingGapInterface:=function(arg)
local
    IT,
    T1, T2, T3, C2, C2p, M, P, C3,
    PGGen, PGGen33, PGMat33, PGMatinv, PGind,
    o0,o1,o2,o3,o4,o5,
    G,Gp,R,CB,
    MatToPow,GapToPow,Fbarhomotopyindinv,
    Homotopydeg1,Homotopydeg2,Homotopydeg3,Homotopydeg4,
    func, funcs,receive,
    Gen1, Gen2, Gen3, Gen4, GensGAP,
    i,j,k,p,x,y;

#####################################################################
MatToPow:=function(mat)            #given 4x4 matrix, output power
local i, mat33, trans;

mat33:=List([1..3],i->List([1..3],j->mat[i,j]));

i:=Position(PGMat33,mat33);

trans:=mat*PGMatinv[i];

return Concatenation(List([1..3],x->trans[x,4]),PGind[i]);
end;

#####################################################################
GapToPow:=function(i)            #given the index i s.t. mat:=R!.elts[i], output the power

return MatToPow(TransposedMat(PreImage(Gp,R!.elts[i])));
end;
#####################################################################
Fbarhomotopyindinv:=function(i,lst)            #This is the function Fbarhomotopyindinv in Mathematica

return List(lst,x->Concatenation(x,[i]));
end;
#####################################################################


    ####################BEGIN TO READ THE INPUT##################

IT:=arg[1];

T1:=[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]]; #standard translation T1
T2:=[[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]]; #standard translation T2
T3:=[[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]]; #standard translation T3


PGMat33 := [];
PGMatinv := [];
PGind := [];

#below are to be read from the file:
Read("~/Downloads/Space_Group_Cocycles.gi");


if IT=47 then
    PGGen:=[M147,M247,M347];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in47,Am2in47,Am3in47,Axin47,Ayin47,Azin47],[],[]];
elif IT=48 then
    PGGen:=[M148,M248,M348];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in48,Am2in48,Am3in48,Axyzin48],[],[]];
elif IT=49 then
    PGGen:=[M149,M249,M349];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in49,Am2in49,Am3in49,Axin49,Ayin49],[],[]];
elif IT=50 then
    PGGen:=[M150,M250,M350];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in50,Am2in50,Am3in50,Azin50],[BGAPin50],[]];
elif IT=51 then
    PGGen:=[M151,M251,M351];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in51,Am2in51,Am3in51,Ayin51,Azin51],[],[]];
elif IT=52 then
    PGGen:=[M152,M252,M352];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in52,Am2in52,Am3in52],[B1in52,B2in52],[]];
elif IT=53 then
    PGGen:=[M153,M253,M353];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in53,Am2in53,Am3in53,Ayin53],[Bin53],[]];
elif IT=54 then
    PGGen:=[M154,M254,M354];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in54,Am2in54,Am3in54,Ayin54],[],[]];
elif IT=55 then
    PGGen:=[M155,M255,M355];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in55,Am2in55,Am3in55,Azin55],[Bin55],[]];
elif IT=56 then
    PGGen:=[M156,M256,M356];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in56,Am2in56,Am3in56],[B1in56,B2in56],[]];
elif IT=57 then
    PGGen:=[M157,M257,M357];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in57,Am2in57,Am3in57,Axin57],[],[]];
elif IT=58 then
    PGGen:=[M158,M258,M358];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in58,Am2in58,Am3in58],[B1in58,B2in58,B3in58],[CGAPin58]];
elif IT=59 then
    PGGen:=[M159,M259,M359];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in59,Am2in59,Am3in59,Azin59],[Bin59],[]];
elif IT=60 then
    PGGen:=[M160,M260,M360];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in60,Am2in60,Am3in60],[Bin60],[]];
elif IT=61 then
    PGGen:=[M161,M261,M361];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in61,Am2in61,Am3in61],[],[CGAPin61]];
elif IT=62 then
    PGGen:=[M162,M262,M362];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in62,Am2in62,Am3in62],[B1in62,B2in62],[]];
elif IT=63 then
    PGGen:=[M163,M263,M363];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in63,Am2in63,Am3in63,Axyin63],[Bin63],[]];
elif IT=64 then
    PGGen:=[M164,M264,M364];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in64,Am2in64,Am3in64,Axyin64],[],[CGAPin64]];
elif IT=65 then
    PGGen:=[M165,M265,M365];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in65,Am2in65,Am3in65,Axyin65,Azin65],[Bin65],[]];
elif IT=66 then
    PGGen:=[M166,M266,M366];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in66,Am2in66,Am3in66,Axyin66],[B1in66,B2in66],[]];
elif IT=67 then
    PGGen:=[M167,M267,M367];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in67,Am2in67,Am3in67,Axyin67,Azin67],[],[]];
elif IT=68 then
    PGGen:=[M168,M268,M368];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in68,Am2in68,Am3in68,Axyin68],[],[CGAPin68]];
elif IT=69 then
    PGGen:=[M169,M269,M369];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in69,Am2in69,Am3in69,Axyin69,Axzin69],[],[CGAPin69]];
elif IT=70 then
    PGGen:=[M170,M270,M370];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in70,Am2in70,Am3in70],[Bin70],[CGAPin70]];
elif IT=71 then
    PGGen:=[M171,M271,M371];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in71,Am2in71,Am3in71,Axyzin71],[B1in71,B2in71,B3in71],[CGAPin71]];
elif IT=72 then
    PGGen:=[M172,M272,M372];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Am1in72,Am2in72,Am3in72,Axyzin72],[Bin72],[]];
elif IT=73 then
    PGGen:=[M173,M273,M373];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in73,Am2in73,Am3in73,Axyzin73],[],[]];
elif IT=74 then
    PGGen:=[M174,M274,M374];
    G:=Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Am1in74,Am2in74,Am3in74,Axyzin74],[B1in74,B2in74],[]];
elif IT=91 then
    PGGen:=[C291, C2p91, C2pp91];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin91,Acppin91,Axyin91],[BGAPin91],[]];
elif IT=92 then
    PGGen:=[C292, C2p92, C2pp92];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin92,Acppin92],[BGAPin92],[CGAPin92]];
elif IT=93 then
    PGGen:=[C293, C2p93, C2pp93];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin93,Acppin93,Axyin93,Azin93],[BGAPin93],[]];
elif IT=94 then
    PGGen:=[C294, C2p94, C2pp94];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin94,Acppin94,Azin94],[BGAPin94],[CGAPin94]];
elif IT=95 then
    PGGen:=[C295, C2p95, C2pp95];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin91,Acppin91,Axyin91],[BGAPin91],[]];
#elif IT=96 then
#    PGGen:=[C296, C2p96, C2pp96]; #incompatible with the cocycles of No.92!!
#    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
#    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
#    funcs:=[[Acpin92,Acppin92],[BGAPin92],[CGAPin92]];
elif IT=97 then
    PGGen:=[C297, C2p97, C2pp97];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin97,Acppin97,Axyzin97],[Bdeltain97,B2in97,B3in97],[C1in97,C2in97]];
elif IT=98 then
    PGGen:=[C298, C2p98, C2pp98];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin98,Acppin98,Axyzin98],[BGAPin98],[CGAPin98]];
elif IT=102 then
    PGGen:=[C2102, Mp102, M102];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Ampin102,Amin102,Axyzin102],[Bdeltain102],[CGAPin102]];
elif IT=106 then
    PGGen:=[C2106, Mp106, M106];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Ampin106,Amin106],[Bdeltain106,B2in106,B3in106],[C1in106,C2in106]];
elif IT=109 then
    PGGen:=[C2109, Mp109, M109];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Ampin109,Amin109],[Bdeltain109,B2in109],[C1in109,C2in109]];
elif IT=110 then
    PGGen:=[C2110, Mp110, M110];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Ampin110,Amin110],[Bdeltain110],[CGAPin110]];
elif IT=111 then
    PGGen:=[C2111, C2p111, M111];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin111,Amin111,Axyin111,Azin111],[Bdeltain111,Bxyin111],[]];
elif IT=112 then
    PGGen:=[C2112, C2p112, M112];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin112,Amin112,Axyin112],[Bdeltain112,Bxyin112,Bczin112,Bzxyin112],[]];
elif IT=113 then
    PGGen:=[C2113, C2p113, M113];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin113,Amin113,Azin113],[Bdeltain113,Bcxyin113],[Cin113]];
elif IT=114 then
    PGGen:=[C2114, C2p114, M114];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin114,Amin114],[Bdeltain114,Bcxyin114,Bcxyzin114],[C1in114,C2in114]];
elif IT=115 then
    PGGen:=[C2115, C2p115, M115];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin115,Amin115,Axyin115,Azin115],[Bdeltain115,Bxyin115],[]];
elif IT=116 then
    PGGen:=[C2116, C2p116, M116];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin116,Amin116,Axyin116],[Bdeltain116,Bxyin116,Bcpzin116],[]];
elif IT=117 then
    PGGen:=[C2117, C2p117, M117];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin117,Amin117,Azin117],[Bdeltain117,B2in117],[CGAPin117]];
elif IT=118 then
    PGGen:=[C2118, C2p118, M118];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Acpin118,Amin118,Axyzin118],[Bdeltain118],[CGAPin118]];
elif IT=119 then
    PGGen:=[C2119, C2p119, M119];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin119,Amin119,Axyzin119],[Bdeltain119,B2in119,B3in119],[C1in119,C2in119]];
elif IT=120 then
    PGGen:=[C2120, C2p120, M120];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin120,Amin120,Axyzin120],[Bdeltain120,B2in120],[]];
elif IT=121 then
    PGGen:=[C2121, C2p121, M121];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Acpin121,Amin121,Axyzin121],[Bdeltain121,Bphiin121,Bxyzin121],[CGAPin121]];
elif IT=123 then
    PGGen:=[C2123, C2p123, M123, P123];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin123,Amin123,Acpin123,Axyin123,Azin123],[Bdeltain123,Bxyin123],[]];
elif IT=124 then
    PGGen:=[C2124, C2p124, M124, P124];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin124,Amin124,Acpin124,Axyin124],[Bdeltain124,Bxyin124],[]];
elif IT=125 then
    PGGen:=[C2125, C2p125, M125, P125];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin125,Amin125,Acpin125,Azin125],[Bdeltain125,Bcxyin125],[]];
elif IT=127 then
    PGGen:=[C2127, C2p127, M127, P127];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin127,Amin127,Acpin127,Azin127],[Bdeltain127,Bcxyin127],[CGAPin127]];
elif IT=128 then
    PGGen:=[C2128, C2p128, M128, P128];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin128,Amin128,Acpin128],[Bdeltain128,Bcxyin128,Bcxyzin128],[C1GAPin128,C2GAPin128]];
elif IT=129 then
    PGGen:=[C2129, C2p129, M129, P129];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin129,Amin129,Acpin129,Azin129],[Bdeltain129,Bcxyin129],[]];
elif IT=130 then
    PGGen:=[C2130, C2p130, M130, P130];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin130,Amin130,Acpin130],[Bdeltain130,Bcxyin130],[]];
elif IT=131 then
    PGGen:=[C2131, C2p131, M131, P131];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin131,Amin131,Acpin131,Axyin131],[Bdeltain131,Bczin131,Bxyin131,Bzxyin131],[]];
elif IT=132 then
    PGGen:=[C2132, C2p132, M132, P132];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin132,Amin132,Acpin132,Axyin132],[Bdeltain132,Bmzin132,Bxyin132],[]];
elif IT=134 then
    PGGen:=[C2134, C2p134, M134, P134];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin134,Amin134,Acpin134,Axyzin134],[Bdeltain134],[]];
elif IT=135 then
    PGGen:=[C2135, C2p135, M135, P135];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin135,Amin135,Acpin135],[Bdeltain135,Bcxyin135,Bczin135],[CGAPin135]];
elif IT=136 then
    PGGen:=[C2136, C2p136, M136, P136];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin136,Amin136,Acpin136],[Bdeltain136,Bcxyin136,Bpxyzin136,Bmzin136],[CGAP1in136,CGAP2in136]];
elif IT=138 then
    PGGen:=[C2138, C2p138, M138, P138];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin138,Amin138,Acpin138],[Bdeltain138,Bmzin138,Bcxyin138,Bpxyzin138],[]];
elif IT=139 then
    PGGen:=[C2139, C2p139, M139, P139];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Aiin139,Amin139,Acpin139,Axyzin139],[Bdeltain139,Bphiin139,Bxyzin139],[CGAPin139]];
elif IT=140 then
    PGGen:=[C2140, C2p140, M140, P140];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));
    funcs:=[[Aiin140,Amin140,Acpin140,Axyzin140],[Bdeltain140,Bcxyin140],[]];
elif IT=142 then
    PGGen:=[C2142, C2p142, M142, P142];
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4]);
    Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
    funcs:=[[Aiin142,Amin142,Acpin142],[Bdeltain142],[]];
fi;


o0:=Length(PGGen);



if o0 = 1 then
    G:= Group(T1,T2,T3,PGGen[1]);
    PGGen33 := [List([1..3],i->List([1..3],j->PGGen[1][i,j]))];
    for o1 in [0..(Order(PGGen33[1])-1)] do
        Append(PGind,[[o1]]);
        Append(PGMat33,[PGGen33[1]^o1]);
        Append(PGMatinv,[(PGGen[1]^o1)^(-1)]);
    od;
fi;
if o0 = 2 then
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2]);
    PGGen33 := List([1..2],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..(Order(PGGen33[2])-1)] do
            Append(PGind,[[o1,o2]]);
            Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2]);
            Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2)^(-1)]);
        od;
    od;
fi;
if o0 = 3 then
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3]);
    PGGen33 := List([1..3],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..(Order(PGGen33[2])-1)] do
            for o3 in [0..(Order(PGGen33[3])-1)] do
                Append(PGind,[[o1,o2,o3]]);
                Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3]);
                Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3)^(-1)]);
            od;
        od;
    od;
fi;
if o0 = 4 then
    PGGen33 := List([1..4],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..(Order(PGGen33[2])-1)] do
            for o3 in [0..(Order(PGGen33[3])-1)] do
                for o4 in [0..(Order(PGGen33[4])-1)] do
                    Append(PGind,[[o1,o2,o3,o4]]);
                    Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3*PGGen33[4]^o4]);
                    Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3*PGGen[4]^o4)^(-1)]);
                od;
            od;
        od;
    od;
fi;
if o0 = 5 then
    G:= Group(T1,T2,T3,PGGen[1],PGGen[2],PGGen[3],PGGen[4],PGGen[5]);
    PGGen33 := List([1..5],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..(Order(PGGen33[2])-1)] do
            for o3 in [0..(Order(PGGen33[3])-1)] do
                for o4 in [0..(Order(PGGen33[4])-1)] do
                    for o5 in [0..(Order(PGGen33[5])-1)] do
                        Append(PGind,[[o1,o2,o3,o4,o5]]);
                        Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3*PGGen33[4]^o4*PGGen33[5]^o5]);
                        Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3*PGGen[4]^o4*PGGen[5]^o5)^(-1)]);
                    od;
                od;
            od;
        od;
    od;
fi;




R:=ResolutionAlmostCrystalGroup(Image(Gp),7);

Homotopydeg1:=List([1..R!.dimension(1)],x->List(R!.boundary(1,x),y->[y[2]]));
Homotopydeg2:=List([1..R!.dimension(2)],x->Concatenation(List(R!.boundary(2,x),y->Fbarhomotopyindinv(y[2],Homotopydeg1[AbsInt(y[1])]))));
Homotopydeg3:=List([1..R!.dimension(3)],x->Concatenation(List(R!.boundary(3,x),y->Fbarhomotopyindinv(y[2],Homotopydeg2[AbsInt(y[1])]))));


CB:=[];
for p in [1..3] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;


Gen1:=[];
for func in funcs[1] do
    Append(Gen1,[CB[1].cocycleToClass(List([1..R!.dimension(1)],x->RemInt(Sum(List(Homotopydeg1[x],y->func(GapToPow(y[1])))),2)))]);
od;

Gen2:=[];
for func in funcs[2] do
    Append(Gen2,[CB[2].cocycleToClass(List([1..R!.dimension(2)],x->RemInt(Sum(List(Homotopydeg2[x],y->func(GapToPow(y[1]),GapToPow(y[2])))),2)))]);
od;

Gen3:=[];
for func in funcs[3] do
    Append(Gen3,[CB[3].cocycleToClass(List([1..R!.dimension(3)],x->RemInt(Sum(List(Homotopydeg3[x],y->func(GapToPow(y[1]),GapToPow(y[2]),GapToPow(y[3])))),2)))]);
od;

Gen4:=[];


GensGAP:=Mod2RingGenerators(R,6,3);


Print("Degree-1 generator:", Gen1,"=",GensGAP[1],"\n");
Print("Degree-2 generator:", Gen2,"=",GensGAP[2],"\n");
Print("Degree-3 generator:", Gen3,"=",GensGAP[3],"\n");

if (IT in [108, 109, 120, 130, 136, 140, 142, 197, 204, 230]) = true then
    Gen4 := GensGAP[4];
fi;

Print("Degree-4 generator:", GensGAP[4],"\n");
Print("Degree-5 generator:", GensGAP[5],"\n");
Print("Degree-6 generator:", GensGAP[6],"\n");


Mod2RingGensAndRels(IT,3,R,[Gen1,Gen2,Gen3,Gen4]);



return true;
end;
#####################################################################
#####################################################################
