
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
        Gens, GensLett, Cupped, CupRelsLett,
        CupBase2all, CupBase2, CupBase3,CupBase4,CupBase5,CupBase6,
        CupBase2Lett,CupBase3Lett,CupBase4Lett,CupBase5Lett,CupBase6Lett,
        CupRel2Lett, CupRel3Lett, CupRel4Lett, CupRel5Lett, CupRel6Lett,
        cupped,
        Lett1, Lett2, IO,
        Letters, AddLetters,
        uCocycle, vCocycle, uvCocycle, uChainMap, ww,
        sol, solrel, cc, CB, CohomologyBasis, TR,
        BasisP, BasisQ, SmithRecord, GF2ToZ, IToPosition,
        NonNegativeVec, Letter2Monomial, PrintMonomialString,
        i, p, q, r, u, v, ln, iu, iv, w, x, sw;

#Standard input: arg[1] = IT (# of space group), arg[2] = n (relations up to deg(n) is calculated)
#e.g.: Mod2RingGensAndRels(89);
#e.g.: Mod2RingGensAndRels(89,6);




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

Print("===========================================\n");
Print("Begin Group No. ", IT, ":\n");

if IT = 219 then
    Print("Caution: the group being calculated (No. 219) has two degree-6 generators!\n");
    Print("These degree-6 generators have not been included in the current code, ");
    Print("the higher degree (degree-12) relations have not been worked out.\n");
fi;
if IT = 226 or IT = 228 then
    Print("Caution: the group being calculated (either No. 226 or No. 228) has a degree-6 generator!\n");
    Print("This degree-6 generator has not been included in the current code, ");
    Print("as the length-7 resolution required to get the degree-6 generator exceeds the current memory limit.\n");
    Print("the higher degree (degree-12) relations have not been worked out.\n");
fi;

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



Letters := [["A1","A2","A3","A4","A5","A6","A7"],["B1","B2","B3","B4","B5","B6","B7"],                             ["C1","C2","C3","C4","C5","C6","C7"],["D1","D2","D3","D4","D5","D6","D7"]];


Gens:=Mod2RingGenerators(R,4,spacedim);
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

GensLett:=CohomologyBasis(List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))],i->1));
#later we will want to use any generators that we designate.
#GensLett records the powers, later to be used to convert generators/relations to letters.

CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

Print("Number of Generators at degrees 1-4: ");
Print([Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4)]);
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

iu :=1;
for u in Gen1 do
iv :=1;
for v in Gen1 do
    Lett1 := GensLett[iu]+GensLett[iv];
    if iv>=iu then
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
                if sol = fail then                                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase2,[cupped]);
                    Append(CupBase2Lett,[Lett1]);
                else                              #if u-cup-v is expressable by other cocycles
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

#Begins: degree-2 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase2 do
iv :=1;
for v in Gen1 do
    Lett1 := CupBase2Lett[iu] + GensLett[iv];
    IO := false;
    for Lett2 in CupRelsLett do
        if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
        IO :=true;       #then rewrite IO
        fi;
    od;
    if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        cupped :=Mod2CupProduct(R,u,v,2,1,CB[2],CB[1],CB[3]);      #then calculate u-cup-v
        
        if cupped = List([1..Cohomology(TR,3)],x->0) then         #if u-cup-v is a coboundary
            Append(CupRelsLett,[Lett1]);
            Append(CupRel3Lett,[[Lett1]]);
        else
            if CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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

#Step-1 begins here: degree-3 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase3 do
iv :=1;
for v in Gen1 do
    Lett1 := CupBase3Lett[iu] + GensLett[iv];
    IO := false;
    for Lett2 in CupRelsLett do
        if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
        IO :=true;       #then rewrite IO
        fi;
    od;
    if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        cupped :=Mod2CupProduct(R,u,v,3,1,CB[3],CB[1],CB[4]);      #then calculate u-cup-v
        
        if cupped = List([1..Cohomology(TR,4)],x->0) then         #if u-cup-v is a coboundary
            Append(CupRelsLett,[Lett1]);
            Append(CupRel4Lett,[[Lett1]]);
        else
            if CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
    fi;
iv := iv+1;
od;
iu := iu+1;
od;
#
#
#Step-1 finished: degree-3 cup with degree-1-gen

#Step-2 begins here: degree-2-gen cup with degree-2-gen
#
#
iu :=1;
for u in Gen2 do
iv :=1;
for v in Gen2 do
    if iv>=iu then
        Lett1 := GensLett[Length(Gen1)+iu] + GensLett[Length(Gen1)+iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            IO :=true;       #then rewrite IO
            fi;
        od;
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
            cupped :=Mod2CupProduct(R,u,v,2,2,CB[2],CB[2],CB[4]);      #then calculate u-cup-v
        
            if cupped = List([1..Cohomology(TR,4)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel4Lett,[[Lett1]]);
            else
                if CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
        Lett1 := CupBase4Lett[iu] + GensLett[iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            IO :=true;       #then rewrite IO
            fi;
        od;
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
    
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
        
            if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel5Lett,[[Lett1]]);
            else
                if CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
        Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            IO :=true;       #then rewrite IO
            fi;
        od;
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        
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
        
            if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel5Lett,[[Lett1]]);
            else
                if CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
        Lett1 := CupBase5Lett[iu] + GensLett[iv];
        IO := false;
        for Lett2 in CupRelsLett do
            if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            IO :=true;       #then rewrite IO
            fi;
        od;
        if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
    
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
        
            if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            else
                if CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];
            IO := false;
            for Lett2 in CupRelsLett do
                if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;       #then rewrite IO
                fi;
            od;
            if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
    
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
        
                if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel6Lett,[[Lett1]]);
                else
                    if CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
            fi;
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen


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
            IO := false;
            for Lett2 in CupRelsLett do
                if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
                IO :=true;       #then rewrite IO
                fi;
            od;
            if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        
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
                else
                    if CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
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
