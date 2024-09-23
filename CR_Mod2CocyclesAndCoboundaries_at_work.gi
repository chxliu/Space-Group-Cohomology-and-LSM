LoadPackage("HAP");

#Read("~/Downloads/Space_Group_Cocycles.gi");


#Functions in this module:
#
#GF2ToZ
#Letter2Monomial
#PrintMonomialString
#CR_Mod2CocyclesAndCoboundaries
#Mod2CupProduct
#Mod2RingGenerators
#Mod2RingGensAndRels
#PointGroupTranslationExtension
#IrreducibleWyckoffPoints
#SpaceGroupCohomologyRingGapInterface



#####################################################################
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
######################################################################


#####################################################################
#####################################################################

Letter2Monomial:=function(vec,GensDim,name)
local j,k,u,v,lett;

lett := [];
u := 1;
if name = [] then
    for j in [1..Length(GensDim)] do
        if GensDim[j] = 1 then
            if vec[u]>1 then
                Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],"^",vec[u]],"")]);
            elif vec[u]=1 then
                Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j]],"")]);
            fi;
            u := u+1;
        else
            for k in [1..GensDim[j]] do
                if vec[u]>1 then
                    Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k),"^",vec[u]],"")]);
                elif vec[u]=1 then
                    Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k)],"")]);
                fi;
                u := u+1;
            od;
        fi;
    od;
else
    for j in [1..Length(GensDim)] do
    
        for k in [1..GensDim[j]] do
            if vec[u]>1 then
                Append(lett,[JoinStringsWithSeparator([name[u],"^",vec[u]],"")]);
            elif vec[u]=1 then
                Append(lett,[JoinStringsWithSeparator([name[u]],"")]);
            fi;
            u := u+1;
        od;
    od;
fi;
return JoinStringsWithSeparator(lett, ".");
end;

#####################################################################
######################################################################


#####################################################################
#####################################################################

PrintMonomialString:=function(arg)
local poly,vecs,GensDim,sep,gennames,space;

vecs := arg[1];
GensDim := arg[2];
sep := arg[3];
gennames := [];
if Length(arg) >=4 then
    gennames := arg[4];
fi;
space:="  ";
if Length(arg) >=5 then
    space := arg[5];
fi;

poly:=List(vecs,x->Letter2Monomial(x,GensDim,gennames));
Print(JoinStringsWithSeparator(poly,sep),space);
return 0;
end;

######################################################################
#####################################################################


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
        BasisP, BasisQ,
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
        if Cupped = [] then
            Gens := Cups*Z(2);
        else
            CuppedBasis := List(BaseMat(Cupped),ShallowCopy);
            Gens := BaseSteinitzVectors(Cups*Z(2),CuppedBasis)!.factorspace;
        fi;

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
        R,n,GG,IT,Gen0,Gen1,Gen2,Gen3,Gen4,Gen5,Gen6,ss6b,ss6,spacedim,GenDim1to4,GenDeg1to4,
        Gens, GensLett, Cupped, CupRelsLett, CupTemp, CupTempLett,
        CupBase2all,CupBase2,CupBase3,CupBase4,CupBase5,CupBase6,CupBase7,CupBase8,
        CupBase1Lett, CupBase2Lett,CupBase3Lett,CupBase4Lett,CupBase5Lett,CupBase6Lett,CupBase7Lett,CupBase8Lett,
        CupRel2Lett, CupRel3Lett, CupRel4Lett, CupRel5Lett, CupRel6Lett, CupRel7Lett, CupRel8Lett,
        RelReduceLett, RelReduceMat, RelReduceVec, RelRedLen, solrelred,
        cupped,M1,sum,row,BasisImaged2,cuppedRaw,CupBase6RawCobandCoc,CupBase6Raw,
        Lett1, Lett2, mono, IO,
        uCocycle, vCocycle, uvCocycle, uChainMap, ww,
        sol, sol1, solrel, cc, CB, CohomologyBasis, TR,
        BasisP, BasisQ, SmithRecord, IToPosition,
        #NonNegativeVec,
        i,j,p,q,r,s,t,u,v,w,x,y,z, ln, rk, rk1, ip,iq,ir,is,it,iu,iv,iw,ix,iy,iz,sw;

#Standard input: arg[1] = IT (# of space group), arg[2] = n (relations up to deg(n) is calculated)
#e.g.: Mod2RingGensAndRels(89);
#e.g.: Mod2RingGensAndRels(89,3);
#e.g.: Mod2RingGensAndRels(89,3,R89);
#e.g.: Mod2RingGensAndRels(89,3,R89,Gens);




if Length(arg)=1 then
    spacedim:=3;
fi;

if Length(arg)>=2 then
    spacedim:=arg[2];
fi;


IT := arg[1];  #Group number in International Table for Crystallography (ITC)

Print("===========================================\n");
#Print("Begin Group No. ", IT, ":\n");



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
Gen5:=[]; #based on the posteriori fact that no space group has degree-5 generators.
Gen6:=[];

if IT = 219 then
    #Print("Caution: the group being calculated (No. 219) has two degree-6 generators!\n");
    ss6b :=CR_Mod2CocyclesAndCoboundaries(R,6)!.cocyclesBasis;
    ss6 :=CR_Mod2CocyclesAndCoboundaries(R,6);
    #Print([ss6!.cocycleToClass(ss6b[Length(ss6b)-1]),ss6!.cocycleToClass(ss6b[Length(ss6b)])]);
    Append(Gen6,[ss6!.cocycleToClass(ss6b[Length(ss6b)-1]),ss6!.cocycleToClass(ss6b[Length(ss6b)])]);
fi;
if IT = 226 or IT=228 then
    Print("Caution: the group being calculated (Either No. 226 or No. 228) has degree-6 generator(s)!\n");
    Print("They cannot be obtained in the current program, as the length-7 resolution required to get the degree-6 generator exceeds the current memory limit.\n");
    #ss6b :=CR_Mod2CocyclesAndCoboundaries(R,6)!.cocyclesBasis;
    #ss6 :=CR_Mod2CocyclesAndCoboundaries(R,6);
    #Append(Gen6,[ss6!.cocycleToClass(ss6b[Length(ss6b)])]);
fi;

GenDim1to4:=[Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4)];
GenDeg1to4:=Concatenation(List([1..Length(Gen1)],x->1),List([1..Length(Gen2)],x->2),List([1..Length(Gen3)],x->3),List([1..Length(Gen4)],x->4));

if Length(Gen4)>0 or Length(Gen6)>0 then
    #Print("Caution: this group contains degree-4 generators!\n");
    Print("Containing higher (>6) degree relations -- needs to be calculated manually!! \n");
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
#NonNegativeVec:=function(v)
#local ps,k;
#
#ps := true;
#for k in [1..Length(v)] do
#    if v[k] < 0 then
#    ps := false;
#    fi;
#od;
#return ps;
#end;
######################################################################


GensLett:=CohomologyBasis(List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))],i->1));
#GensLett records the powers, later to be used to convert generators/relations to letters.

Gen0 := List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))],i->0);

TR:=HomToIntegersModP(R,2);            #Apply the Hom functor.

#if Length(Gen4) = 0 then
#    n:=6; #This is the highest degree at which the relations are calculated. For now we fix it to be 6.
#else
#    n:=8;
#fi;

n:=Length(Size(R))-1;     #This is the highest degree at which the relations are calculated for groups No. <= 220.
                          #for 221, Length(Size(R)) = 6, n = 5, and we will use the "raw code" to calculate relations at degree 6.


CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

#Print("Number of Generators at degrees 1-4: ");
#Print(GenDim1to4);

####################### r = 1 ##########################

#Print("Chosen basis at degree 1:\n");
#PrintMonomialString(List([1..Length(Gen1)],x->GensLett[x]),GenDim1to4,",");

#Print("Matching dimensions... [dim(Chosen basis) = dim(H)]\n");
#Print("dim(H^1)=", Cohomology(TR,1),", ");

CupBase1Lett:= List([1..Length(Gen1)],x->GensLett[x]);


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
                    else                       #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                        if (Lett1 in CupBase2Lett) = false then
                            Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                            Append(CupRel2Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
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
                else                              #if u-cup-u is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                    if (Lett1 in CupBase2Lett) = false then
                        Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                        Append(CupRel2Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    fi;
                fi;
            fi;
        fi;
    iu := iu+1;
od;

#Append(CupBase2,Gen2);
for iu in [(1+Length(Gen1))..(Length(Gen1)+Length(Gen2))] do
    if CupBase2 = [] then
        Append(CupBase2,[Gen2[iu-Length(Gen1)]]);
        Append(CupBase2Lett,[GensLett[iu]]);
    else
        sol :=SolutionMat(CupBase2*Z(2),Gen2[iu-Length(Gen1)]*Z(2));
        if sol = fail then
            Append(CupBase2,[Gen2[iu-Length(Gen1)]]);
            Append(CupBase2Lett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-2 generator(s)!!\n");
        fi;
    fi;
    #Append(CupBase2Lett,[GensLett[iu]]);
od;


#Print("Independent relations at degree 2:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 2:\n");
#PrintMonomialString(CupBase2Lett,GenDim1to4,",");

if Length(CupBase2Lett) = Cohomology(TR,2) then
    Print("");#Print("dim(H^2)=", Cohomology(TR,2),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^2) = ", Length(CupBase2Lett) - Cohomology(TR,2),"\n");
fi;

#Print("with relations:\n");
#Print(CupRelsLett,"\n");
##Begin printing the relations at deg 2:
if Length(CupRel2Lett) > 0 then
    #Print(Length(CupRel2Lett)," at deg 2: ");List(CupRel2Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R2:  ");List(CupRel2Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
##End printing the relations at deg 2.


####################### r = 3 ##########################


CupBase3 :=[];
CupBase3Lett :=[];
CupRel3Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];



#### Begin preparation for relation reduction
####
iu := 1;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            if (p+q+r)*GenDeg1to4 = 3 then
                Append(RelReduceLett,[p+q+r]);
                iu := iu+1;
            fi;
        od;
    od;
od;

RelRedLen := iu;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix

for u in CupBase1Lett do
    for v in CupRel2Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        Append(RelReduceMat,[RelReduceVec]);
    od;
od;
####
#### End preparation for relation reduction


#Begins: degree-2 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase2 do
    iv :=1;
    for v in Gen1 do
        
        cupped :=Mod2CupProduct(R,u,v,2,1,CB[2],CB[1],CB[3]);      #then calculate u-cup-v
        Lett1 := CupBase2Lett[iu] + GensLett[iv];
                
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        #Lett1 := CupBase2Lett[iu] + GensLett[iv];
        #IO := false;
        #for Lett2 in CupRelsLett do
        #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
        #        IO :=true;                      #then rewrite IO
        #        Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
        #        Append(CupTempLett,[Lett1]);
        #    fi;
        #od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations

        #if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,


        #### Start: determine whether u-cup-v goes to the basis
        ####
        if (cupped = List([1..Cohomology(TR,3)],x->0)) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
        else
            if CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase3,[cupped]);
                Append(CupBase3Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase3*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase3,[cupped]);
                    Append(CupBase3Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase3Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
                    
                    
        #Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
        #Append(CupTempLett,[Lett1]);
                
                
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase3Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
            
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel3Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Finished: degree-2 cup with degree-1-gen


#Append(CupBase3,Gen3);
for iu in [(1+Length(Gen1)+Length(Gen2))..(Length(Gen1)+Length(Gen2)+Length(Gen3))] do
    if CupBase3 = [] then
        Append(CupBase3,[Gen3[iu-(Length(Gen1)+Length(Gen2))]]);
        Append(CupBase3Lett,[GensLett[iu]]);
    else
        sol :=SolutionMat(CupBase3*Z(2),Gen3[iu-(Length(Gen1)+Length(Gen2))]*Z(2));
        if sol = fail then
            Append(CupBase3,[Gen3[iu-(Length(Gen1)+Length(Gen2))]]);
            Append(CupBase3Lett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-3 generator(s)!!\n");
        fi;
    fi;
    #Append(CupBase3Lett,[GensLett[iu]]);
od;



#if CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase3,[CupTemp[1]]);
#            Append(CupBase3Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;

#if (CupBase3 = []) = false then
#
#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase3*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase3,[cupped]);
#            Append(CupBase3Lett,[Lett1]);
#        fi;
#    od;
#fi;

#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 3: ");
#PrintMonomialString(CupBase3Lett,GenDim1to4,",");
#Print("\n");

if Length(CupBase3Lett) = Cohomology(TR,3) then
    Print("");#Print("dim(H^3)=", Cohomology(TR,3),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^3) = ", Length(CupBase3Lett) - Cohomology(TR,3),"\n");
fi;



##Begin printing the relations at deg 3:
if Length(CupRel3Lett) > 0 then
    #Print(Length(CupRel3Lett)," at deg 3: ");List(CupRel3Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R3:  ");List(CupRel3Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
##End printing the relations at deg 3.




####################### r = 4 ##########################

CupBase4 :=[];
CupBase4Lett :=[];
CupRel4Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];

#### Begin preparation for relation reduction
####
iu := 1;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                if (p+q+r+s)*GenDeg1to4 = 4 then
                    Append(RelReduceLett,[p+q+r+s]);
                    iu := iu+1;
                fi;
            od;
        od;
    od;
od;

RelRedLen := iu;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix

for u in CupBase1Lett do
    for v in CupRel3Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        #Print("RelReduceVec:",RelReduceVec,"\n");
        Append(RelReduceMat,[RelReduceVec]);
        #Print("RelReduceMat:",RelReduceMat,"\n");
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        if (p+q)*GenDeg1to4 = 2 then
            for v in CupRel2Lett do
                RelReduceVec := List([1..RelRedLen],x->0);
                for w in v do
                    RelReduceVec[Position(RelReduceLett,p+q+w)] := 1;
                od;
                Append(RelReduceMat,[RelReduceVec]);
            od;
        fi;
    od;
od;
####
#### End preparation for relation reduction



#Step-1 begins here: degree-3 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase3 do
    iv :=1;
    for v in Gen1 do

        cupped :=Mod2CupProduct(R,u,v,3,1,CB[3],CB[1],CB[4]);      #then calculate u-cup-v
        Lett1 := CupBase3Lett[iu] + GensLett[iv];
        
        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        #Lett1 := CupBase3Lett[iu] + GensLett[iv];
        #IO := false;
        #for Lett2 in CupRelsLett do
        #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
        #        IO :=true;                      #then rewrite IO
        #        Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
        #        Append(CupTempLett,[Lett1]);
        #    fi;
        #od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations
    
        #if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,

        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,4)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase4,[cupped]);
            Append(CupBase4Lett,[Lett1]);
        else
            sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase4,[cupped]);
                Append(CupBase4Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase4Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel4Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-3 cup with degree-1-gen

#Print("after first round (3-cup-1) the dimension of the basis is: ",Length(CupBase4Lett), "\n");

#if CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase4,[CupTemp[1]]);
#            Append(CupBase4Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;



#if (CupBase4 = []) = false then
#
#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase4*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase4,[cupped]);
#            Append(CupBase4Lett,[Lett1]);
#        fi;
#    od;
#fi;



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
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
                    if (Lett1 in CupBase4Lett) = false then
                        Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                        Append(CupRel4Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
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

#Append(CupBase4,Gen4);
for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))] do
    if CupBase4 = [] then
        Append(CupBase4,[Gen4[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3))]]);
        Append(CupBase4Lett,[GensLett[iu]]);
    else
        sol :=SolutionMat(CupBase4*Z(2),Gen4[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3))]*Z(2));
        if sol = fail then
            Append(CupBase4,[Gen4[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3))]]);
            Append(CupBase4Lett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-4 generator(s)!!\n");
        fi;
    fi;
    #Append(CupBase4Lett,[GensLett[iu]]);
od;



#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 4:\n");
#PrintMonomialString(CupBase4Lett,GenDim1to4,",");

if Length(CupBase4Lett) = Cohomology(TR,4) then
    Print("");#Print("dim(H^4)=", Cohomology(TR,4),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^4) = ", Length(CupBase4Lett) - Cohomology(TR,4),"\n");
fi;


##Begin printing the relations at deg 4:
if Length(CupRel4Lett) > 0 then
    #Print(Length(CupRel4Lett)," at deg 4: ");List(CupRel4Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R4:  ");List(CupRel4Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
##End printing the relations at deg 4.



####################### r = 5 ##########################

CupBase5 :=[];
CupBase5Lett :=[];
CupRel5Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];

#### Begin preparation for relation reduction
####
iu := 1;

for ip in [1..(Sum(GenDim1to4)+1)] do
    p := Concatenation([Gen0],GensLett)[ip];
    for iq in [ip..(Sum(GenDim1to4)+1)] do
        q := Concatenation([Gen0],GensLett)[iq];
        for ir in [iq..(Sum(GenDim1to4)+1)] do
            r := Concatenation([Gen0],GensLett)[ir];
            for is in [ir..(Sum(GenDim1to4)+1)] do
                s := Concatenation([Gen0],GensLett)[is];
                for it in [is..(Sum(GenDim1to4)+1)] do
                    t := Concatenation([Gen0],GensLett)[it];
                    if (p+q+r+s+t)*GenDeg1to4 = 5 then
                        Append(RelReduceLett,[p+q+r+s+t]);
                        iu := iu+1;
                    fi;
                od;
            od;
        od;
    od;
od;

RelRedLen := iu;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix

for u in CupBase1Lett do
    for v in CupRel4Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        Append(RelReduceMat,[RelReduceVec]);
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        if (p+q)*GenDeg1to4 = 2 then
            for v in CupRel3Lett do
                RelReduceVec := List([1..RelRedLen],x->0);
                for w in v do
                    RelReduceVec[Position(RelReduceLett,p+q+w)] := 1;
                od;
                Append(RelReduceMat,[RelReduceVec]);
            od;
        fi;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            if (p+q+r)*GenDeg1to4 = 3 then
                for v in CupRel2Lett do
                     RelReduceVec := List([1..RelRedLen],x->0);
                    for w in v do
                        RelReduceVec[Position(RelReduceLett,p+q+r+w)] := 1;
                    od;
                    Append(RelReduceMat,[RelReduceVec]);
                od;
            fi;
        od;
    od;
od;
####
#### End preparation for relation reduction

rk:=RankMatrix(RelReduceMat*Z(2));


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

        Lett1 := CupBase4Lett[iu] + GensLett[iv];


        #### begin checking if u-cup-v contains letters of the lower-degree relations
        ####
        #Lett1 := CupBase4Lett[iu] + GensLett[iv];
        #IO := false;
        
        #for Lett2 in CupRelsLett do
        #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
        #        IO :=true;                      #then rewrite IO
        #        Append(CupTemp,[cupped]);       #and store u-cup-v, whose cocycle-ness to be checked
        #        Append(CupTempLett,[Lett1]);
        #    fi;
        #od;
        ####
        #### finished checking if u-cup-v contains letters of the lower-degree relations

        
        #if IO = false then   #If u-cup-v does not contain patterns in CupRelsLett,
        
        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase5,[cupped]);
            Append(CupBase5Lett,[Lett1]);
        else
            sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase5Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel5Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
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
        
        Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+iv];
                
        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,5)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase5,[cupped]);
            Append(CupBase5Lett,[Lett1]);
        else
            sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase5Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel5Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-2 finished: degree-3-gen cup with degree-2-gen


#Append(CupBase5,Gen5);
for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))] do
    if CupBase5 = [] then
        Append(CupBase5,[Gen5[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))]]);
        Append(CupBase5Lett,[GensLett[iu]]);
    else
        sol :=SolutionMat(CupBase5*Z(2),Gen5[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))]*Z(2));
        if sol = fail then
            Append(CupBase5,[Gen5[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))]]);
            Append(CupBase5Lett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-5 generator(s)!!\n");
        fi;
    fi;
    #Append(CupBase5Lett,[GensLett[iu]]);
od;


#if CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase5,[CupTemp[1]]);
#            Append(CupBase5Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;


#if (CupBase5 = []) = false then
#
#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase5*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase5,[cupped]);
#            Append(CupBase5Lett,[Lett1]);
#        fi;
#    od;
#fi;


#Print("Independent relations:\n");
#Print(CupRelsLett);
#Print("\n");


#Print("Chosen basis at degree 5:\n");
#PrintMonomialString(CupBase5Lett,GenDim1to4,",");


if Length(CupBase5Lett) = Cohomology(TR,5) then
    Print("");#Print("dim(H^5)=", Cohomology(TR,5),", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^5) = ", Length(CupBase5Lett) - Cohomology(TR,5),"\n");
fi;


##Begin printing the relations at deg 5:
if Length(CupRel5Lett) >0 then
    #Print(Length(CupRel5Lett)," at deg 5: ");List(CupRel5Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R5:  ");List(CupRel5Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
##End printing the relations at deg 5.


####################### r = 6 ##########################

CupBase6 :=[];
CupBase6Lett :=[];
CupRel6Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];
M1:=[];
CupBase6Raw:=[];

#Below starts the "Raw Method" for the degree 6 relations for groups 221-230, part1
#
if Length(Size(R)) = 6 then

    #M1 is Dim(R6) x Dim(R5);

    for i in [1..R!.dimension(6)] do
        row:=[];
        for j in [1..R!.dimension(5)] do
            sum:=0;
            for x in R!.boundary(6,i) do
                if AbsoluteValue(x[1])=j then
                    sum := sum + SignInt(x[1]);
                fi;
            od;
        row[j]:= RemInt(sum,2);
        od;
    M1[i]:=row;
    od;
    
    if M1 = [ ] then
        BasisImaged2:=[];
    else
        BasisImaged2:=BaseMat(TransposedMat(M1)*Z(2));
    fi;

    CupBase6RawCobandCoc:=[];
    for i in [1..Length(BasisImaged2)] do
        Append(CupBase6RawCobandCoc,[BasisImaged2[i]]);
    od;

fi;

#
#Above ends the "Raw Method" for the degree 6 relations for groups 221-230, part1



#### Begin preparation for relation reduction
####
iv := 1;
for ip in [1..(Sum(GenDim1to4)+1)] do
    p := Concatenation([Gen0],GensLett)[ip];
    for iq in [ip..(Sum(GenDim1to4)+1)] do
        q := Concatenation([Gen0],GensLett)[iq];
        for ir in [iq..(Sum(GenDim1to4)+1)] do
            r := Concatenation([Gen0],GensLett)[ir];
            for is in [ir..(Sum(GenDim1to4)+1)] do
                s := Concatenation([Gen0],GensLett)[is];
                for it in [is..(Sum(GenDim1to4)+1)] do
                    t := Concatenation([Gen0],GensLett)[it];
                    for iu in [it..(Sum(GenDim1to4)+1)]  do
                        u := Concatenation([Gen0],GensLett)[iu];
                        if (p+q+r+s+t+u)*GenDeg1to4 = 6 then
                            Append(RelReduceLett,[p+q+r+s+t+u]);
                            iv := iv+1;
                        fi;
                    od;
                od;
            od;
        od;
    od;
od;

RelRedLen := iv;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix


for u in CupBase1Lett do
    for v in CupRel5Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        Append(RelReduceMat,[RelReduceVec]);
    od;
od;



for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        if (p+q)*GenDeg1to4 = 2 then
            for v in CupRel4Lett do
                RelReduceVec := List([1..RelRedLen],x->0);
                for w in v do
                    RelReduceVec[Position(RelReduceLett,p+q+w)] := 1;
                od;
                Append(RelReduceMat,[RelReduceVec]);
            od;
        fi;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            if (p+q+r)*GenDeg1to4 = 3 then
                for v in CupRel3Lett do
                    RelReduceVec := List([1..RelRedLen],x->0);
                    for w in v do
                        RelReduceVec[Position(RelReduceLett,p+q+r+w)] := 1;
                    od;
                    Append(RelReduceMat,[RelReduceVec]);
                od;
            fi;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                if (p+q+r+s)*GenDeg1to4 = 4 then
                    for v in CupRel2Lett do
                        RelReduceVec := List([1..RelRedLen],x->0);
                        for w in v do
                            RelReduceVec[Position(RelReduceLett,p+q+r+s+w)] := 1;
                        od;
                        Append(RelReduceMat,[RelReduceVec]);
                    od;
                fi;
            od;
        od;
    od;
od;


####
#### End preparation for relation reduction


rk:=RankMatrix(RelReduceMat*Z(2));


# Below is the good method when resolution at degree 7 has been constructed:
#
#
#
if Length(Size(R)) > 6 then

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

        Lett1 := CupBase5Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase6,[cupped]);
            Append(CupBase6Lett,[Lett1]);
        else
            sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase6Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel6Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
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
        

            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            if cupped = List([1..Cohomology(TR,6)],x->0) then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase6Lett) = false then
               
                #### begin checking whether the relation is reducible from lower degree ones
                ####
                RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
                for w in solrel do
                    #Print(w,"\n");
                    solrelred := Position(RelReduceLett,w);
                            
                    if solrelred = fail then      #if this relation contains new terms then we find a new relation
                        break;
                    else
                        RelReduceVec[solrelred] := 1;
                    fi;
                od;
                if (solrelred = fail) = false then
                    rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                    if (rk = rk1) = false then
                        solrelred :=fail;
                    fi;
                    #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
                fi;
                ####
                #### end checking whether the relation is reducible from lower degree ones

                if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new     relation
                    Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                    Append(CupRel6Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    Append(RelReduceMat,[RelReduceVec]);
                    rk := rk1;
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen




#if CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase6,[CupTemp[1]]);
#            Append(CupBase6Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;


#if (CupBase6 = []) = false then

#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase6*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase6,[cupped]);
#            Append(CupBase6Lett,[Lett1]);
#        fi;
#    od;
#fi;



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
                else                                #if u-cup-v is expressable by other cocycles in basis
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
for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))] do
    if CupBase6 = [] then
        Append(CupBase6,[Gen6[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))]]);
        Append(CupBase6Lett,[GensLett[iu]]);
    else
        sol :=SolutionMat(CupBase6*Z(2),Gen6[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))]*Z(2));
        if sol = fail then
            Append(CupBase6,[Gen6[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5))]]);
            #Append(CupBase6Lett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-6 generator(s)!!\n");
        fi;
    fi;
    #Append(CupBase6Lett,[GensLett[iu]]);
od;



#Print("Chosen basis at degree 6:\n");
#PrintMonomialString(CupBase6Lett,GenDim1to4,",");

if Length(CupBase6) = Cohomology(TR,6) then
    Print("");#Print("dim(H^6)=", Cohomology(TR,6),".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^6) = ", Length(CupBase6) - Cohomology(TR,6),"\n");
fi;

###########################
###########################
###########################
else    ##Below starts the "Raw Method" for the degree 6 relations for groups 221-230, part2
###########################
###########################
###########################

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
        #cupped := CB[6].cocycleToClass(uvCocycle);
        cuppedRaw := uvCocycle;
        ####
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2 ended ####

        Lett1 := CupBase5Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if (SolutionMat(BasisImaged2*Z(2),cuppedRaw*Z(2)) = fail) = false then        #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase6Raw,[cuppedRaw]);
            Append(CupBase6Lett,[Lett1]);
            Append(CupBase6RawCobandCoc,[cuppedRaw]);
        else
            sol :=SolutionMat(CupBase6RawCobandCoc*Z(2),cuppedRaw*Z(2));
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                sol1:=List([(Length(BasisImaged2)+1)..Length(sol)],x->sol[x]);
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase6Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel6Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
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
            #cupped := CB[6].cocycleToClass(uvCocycle);
            cuppedRaw := uvCocycle;
            ####
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2 ended ####
        

            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            if (SolutionMat(BasisImaged2*Z(2),cuppedRaw*Z(2)) = fail) = false then        #if u-cup-v is a coboundary
                solrel := [Lett1];
            
            elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else
                sol :=SolutionMat(CupBase6RawCobandCoc*Z(2),cuppedRaw*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6Raw,[cuppedRaw]);
                    Append(CupBase6Lett,[Lett1]);
                    Append(CupBase6RawCobandCoc,[cuppedRaw]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    sol1:=List([(Length(BasisImaged2)+1)..(Length(sol))],x->sol[x]);
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase6Lett) = false then
               
                #### begin checking whether the relation is reducible from lower degree ones
                ####
                RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
                for w in solrel do
                    #Print(w,"\n");
                    solrelred := Position(RelReduceLett,w);
                            
                    if solrelred = fail then      #if this relation contains new terms then we find a new relation
                        break;
                    else
                        RelReduceVec[solrelred] := 1;
                    fi;
                od;
                if (solrelred = fail) = false then
                    rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                    if (rk = rk1) = false then
                        solrelred :=fail;
                    fi;
                    #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
                fi;
                ####
                #### end checking whether the relation is reducible from lower degree ones

                if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new     relation
                    Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                    Append(CupRel6Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    Append(RelReduceMat,[RelReduceVec]);
                    rk := rk1;
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
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
            #cupped := CB[6].cocycleToClass(uvCocycle);
            cuppedRaw := uvCocycle;
            ####
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2 ended ####
            
            #### Start:  determine whether u-cup-v goes to the basis
            ####
            if (SolutionMat(BasisImaged2*Z(2),cuppedRaw*Z(2)) = fail) = false then        #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            
            elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else
                sol :=SolutionMat(CupBase6RawCobandCoc*Z(2),cuppedRaw*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6Raw,[cuppedRaw]);
                    Append(CupBase6Lett,[Lett1]);
                    Append(CupBase6RawCobandCoc,[cuppedRaw]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    sol1:=List([(Length(BasisImaged2)+1)..(Length(sol))],x->sol[x]);
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
                    if (Lett1 in CupBase6Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel6Lett,[solrel]);
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-3-gen cup with degree-3-gen



#Append(CupBase6,Gen6); #But note that for 221-230, since resolution can only go up to deg 6, we cannot get Gen6.



#Print("Chosen basis at degree 6:\n");
#PrintMonomialString(CupBase6Lett,GenDim1to4,",");

if Length(CupBase6Raw) = [62,11,31,26,45,20,19,6,40,7][IT-220] then
    Print("");#Print("dim(H^6)=", [62,11,31,26,45,20,19,6,40,7][IT-220],".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^6) = ", Length(CupBase6Raw) - [62,11,31,26,45,20,19,6,40,7][IT-220],"\n");
fi;

fi;
#
#
#
#         ##Above ends the "Raw Method" for the degree 6 relations for groups 221-230, part2




##Begin printing the relations at deg 6:
if Length(CupRel6Lett) >0 then
    #Print(Length(CupRel6Lett)," at deg 6: ");List(CupRel6Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R6:  ");List(CupRel6Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
##End printing the relations at deg 6.

####################################################################################
#The part of the code below only runs when relations at r=7 and r=8 are needed, which requires resolution at degree 9 constructed:
#
#
#
#The "if" below is the one that starts working on the r=7 and r=8 relations, which is executed only when Gen4 is non-empty!

if Length(Size(R)) >= 9 then

####################### r = 7 ##########################

CupBase7 :=[];
CupBase7Lett :=[];
CupRel7Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];

#### Begin preparation for relation reduction
####
iw := 1;
for ip in [1..(Sum(GenDim1to4)+1)] do
    p := Concatenation([Gen0],GensLett)[ip];
    for iq in [ip..(Sum(GenDim1to4)+1)] do
        q := Concatenation([Gen0],GensLett)[iq];
        for ir in [iq..(Sum(GenDim1to4)+1)] do
            r := Concatenation([Gen0],GensLett)[ir];
            for is in [ir..(Sum(GenDim1to4)+1)] do
                s := Concatenation([Gen0],GensLett)[is];
                for it in [is..(Sum(GenDim1to4)+1)] do
                    t := Concatenation([Gen0],GensLett)[it];
                    for iu in [it..(Sum(GenDim1to4)+1)]  do
                        u := Concatenation([Gen0],GensLett)[iu];
                        for iv in [iu..(Sum(GenDim1to4)+1)]  do
                            v := Concatenation([Gen0],GensLett)[iv];
                            if (p+q+r+s+t+u+v)*GenDeg1to4 = 7 then
                                Append(RelReduceLett,[p+q+r+s+t+u+v]);
                                iw := iw+1;
                            fi;
                        od;
                    od;
                od;
            od;
        od;
    od;
od;

RelRedLen := iw;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix


for u in CupBase1Lett do
    for v in CupRel6Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        Append(RelReduceMat,[RelReduceVec]);
    od;
od;



for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        if (p+q)*GenDeg1to4 = 2 then
            for v in CupRel5Lett do
                RelReduceVec := List([1..RelRedLen],x->0);
                for w in v do
                    RelReduceVec[Position(RelReduceLett,p+q+w)] := 1;
                od;
                Append(RelReduceMat,[RelReduceVec]);
            od;
        fi;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            if (p+q+r)*GenDeg1to4 = 3 then
                for v in CupRel4Lett do
                    RelReduceVec := List([1..RelRedLen],x->0);
                    for w in v do
                        RelReduceVec[Position(RelReduceLett,p+q+r+w)] := 1;
                    od;
                    Append(RelReduceMat,[RelReduceVec]);
                od;
            fi;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                if (p+q+r+s)*GenDeg1to4 = 4 then
                    for v in CupRel3Lett do
                        RelReduceVec := List([1..RelRedLen],x->0);
                        for w in v do
                            RelReduceVec[Position(RelReduceLett,p+q+r+s+w)] := 1;
                        od;
                        Append(RelReduceMat,[RelReduceVec]);
                    od;
                fi;
            od;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                for t in Concatenation([Gen0],GensLett) do
                    if (p+q+r+s+t)*GenDeg1to4 = 5 then
                        for v in CupRel2Lett do
                            RelReduceVec := List([1..RelRedLen],x->0);
                            for w in v do
                                RelReduceVec[Position(RelReduceLett,p+q+r+s+t+w)] := 1;
                            od;
                            Append(RelReduceMat,[RelReduceVec]);
                        od;
                    fi;
                od;
            od;
        od;
    od;
od;


####
#### End preparation for relation reduction


rk:=RankMatrix(RelReduceMat*Z(2));



#Step-1 begins here: degree-6 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase6 do

    #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 1: ####
    ####
    uCocycle:=CB[6].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,6,1);
    ww:=[];
    for i in [1..(R!.dimension(7))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        #cupped :=Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]);      #then calculate u-cup-v
        
        #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(7))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[7].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 2 ended ####

        Lett1 := CupBase6Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,7)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        else
            if CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase7,[cupped]);
                Append(CupBase7Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase7*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase7Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel7Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-6 cup with degree-1-gen


#Step-2 begins here: (degree-3 cup with degree-2) cup with degree-2-gen
#
#
iu :=1;
for u in CupBase5 do
    if List([1..Length(Gen1)],x->CupBase5Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Bi-cup-Cj
        #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 1: ####
        ####
        uCocycle:=CB[5].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,5,2);
        ww:=[];
        for i in [1..(R!.dimension(7))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
            #cupped :=Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]);      #then calculate u-cup-v
        
            #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(7))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[7].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 2 ended ####
        

            Lett1 := CupBase5Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            if cupped = List([1..Cohomology(TR,7)],x->0) then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            else
                if CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else
                    sol :=SolutionMat(CupBase7*Z(2),cupped*Z(2));
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBase7,[cupped]);
                        Append(CupBase7Lett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase7Lett) = false then
               
                #### begin checking whether the relation is reducible from lower degree ones
                ####
                RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
                for w in solrel do
                    #Print(w,"\n");
                    solrelred := Position(RelReduceLett,w);
                            
                    if solrelred = fail then      #if this relation contains new terms then we find a new relation
                        break;
                    else
                        RelReduceVec[solrelred] := 1;
                    fi;
                od;
                if (solrelred = fail) = false then
                    rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                    if (rk = rk1) = false then
                        solrelred :=fail;
                    fi;
                    #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
                fi;
                ####
                #### end checking whether the relation is reducible from lower degree ones

                if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new     relation
                    Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                    Append(CupRel7Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    Append(RelReduceMat,[RelReduceVec]);
                    rk := rk1;
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-3 cup with degree-2) cup with degree-2-gen




#if CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase7,[CupTemp[1]]);
#            Append(CupBase7Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;


#if (CupBase7 = []) = false then

#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase7*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase7,[cupped]);
#            Append(CupBase7Lett,[Lett1]);
#        fi;
#    od;
#fi;



#Step-3 begins here: degree-4-gen cup with degree-3-gen
#
#
iu :=1;
for u in Gen4 do

    #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,3);
    ww:=[];
    for i in [1..(R!.dimension(7))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7])  -- part 1 ended ####

    iv :=1;
    for v in Gen3 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iu] + GensLett[Length(Gen1)+Length(Gen2)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
            #IO := false;
            #for Lett2 in CupRelsLett do
            #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            #    IO :=true;       #then rewrite IO
            #    fi;
            #od;
        
            #cupped :=Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]);      #then calculate u-cup-v
            
            #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 2: ####
            ####
            vCocycle:=CB[3].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(7))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[7].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 2 ended ####
        
            if cupped = List([1..Cohomology(TR,7)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel7Lett,[[Lett1]]);
            elif CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase7,[cupped]);
                Append(CupBase7Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase7*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                    if (Lett1 in CupBase7Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel7Lett,[solrel]);
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
#Step-3 finished: degree-4-gen cup with degree-3-gen

#Append(CupBase7,Gen7);
#Note that there's no Gen7!!
#for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7))] do
#    if CupBase7 = [] then
#        Append(CupBase7,[Gen7[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))]]);
#        Append(CupBase7Lett,[GensLett[iu]]);
#    else
#        sol :=SolutionMat(CupBase6*Z(2),Gen7[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))]*Z(2));
#        if sol = fail then
#            Append(CupBase7,[Gen7[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))]]);
#            #Append(CupBase7Lett,[GensLett[iu]]);
#        else
#            Print("Error: containing fake degree-7 generator(s)!!\n");
#        fi;
#    fi;
#    #Append(CupBase7Lett,[GensLett[iu]]);
#od;



#Print("Chosen basis at degree 7:\n");
#PrintMonomialString(CupBase7Lett,GenDim1to4,",");

if Length(CupBase7) = Cohomology(TR,7) then
    Print("");#Print("dim(H^7)=", Cohomology(TR,7),".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^7) = ", Length(CupBase7) - Cohomology(TR,7),"\n");
fi;



####################### r = 8 ##########################

CupBase8 :=[];
CupBase8Lett :=[];
CupRel8Lett := [];
CupTemp := [];
CupTempLett := [];
RelReduceLett := [];
RelReduceMat  := [];

#### Begin preparation for relation reduction
####
ix := 1;
for ip in [1..(Sum(GenDim1to4)+1)] do
    p := Concatenation([Gen0],GensLett)[ip];
    for iq in [ip..(Sum(GenDim1to4)+1)] do
        q := Concatenation([Gen0],GensLett)[iq];
        for ir in [iq..(Sum(GenDim1to4)+1)] do
            r := Concatenation([Gen0],GensLett)[ir];
            for is in [ir..(Sum(GenDim1to4)+1)] do
                s := Concatenation([Gen0],GensLett)[is];
                for it in [is..(Sum(GenDim1to4)+1)] do
                    t := Concatenation([Gen0],GensLett)[it];
                    for iu in [it..(Sum(GenDim1to4)+1)]  do
                        u := Concatenation([Gen0],GensLett)[iu];
                        for iv in [iu..(Sum(GenDim1to4)+1)]  do
                            v := Concatenation([Gen0],GensLett)[iv];
                            for iw in [iv..(Sum(GenDim1to4)+1)]  do
                                w := Concatenation([Gen0],GensLett)[iw];
                                if (p+q+r+s+t+u+v+w)*GenDeg1to4 = 8 then
                                    Append(RelReduceLett,[p+q+r+s+t+u+v+w]);
                                    ix := ix+1;
                                fi;
                            od;
                        od;
                    od;
                od;
            od;
        od;
    od;
od;

RelRedLen := ix;

RelReduceMat := [List([1..RelRedLen],x->0)];         #This makes sure that RelReduceMat is not an empty matrix


for u in CupBase1Lett do
    for v in CupRel7Lett do
        RelReduceVec := List([1..RelRedLen],x->0);
        for w in v do
            RelReduceVec[Position(RelReduceLett,u+w)] := 1;
        od;
        Append(RelReduceMat,[RelReduceVec]);
    od;
od;



for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        if (p+q)*GenDeg1to4 = 2 then
            for v in CupRel6Lett do
                RelReduceVec := List([1..RelRedLen],x->0);
                for w in v do
                    RelReduceVec[Position(RelReduceLett,p+q+w)] := 1;
                od;
                Append(RelReduceMat,[RelReduceVec]);
            od;
        fi;
    od;
od;


for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            if (p+q+r)*GenDeg1to4 = 3 then
                for v in CupRel5Lett do
                    RelReduceVec := List([1..RelRedLen],x->0);
                    for w in v do
                        RelReduceVec[Position(RelReduceLett,p+q+r+w)] := 1;
                    od;
                    Append(RelReduceMat,[RelReduceVec]);
                od;
            fi;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                if (p+q+r+s)*GenDeg1to4 = 4 then
                    for v in CupRel4Lett do
                        RelReduceVec := List([1..RelRedLen],x->0);
                        for w in v do
                            RelReduceVec[Position(RelReduceLett,p+q+r+s+w)] := 1;
                        od;
                        Append(RelReduceMat,[RelReduceVec]);
                    od;
                fi;
            od;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                for t in Concatenation([Gen0],GensLett) do
                    if (p+q+r+s+t)*GenDeg1to4 = 5 then
                        for v in CupRel3Lett do
                            RelReduceVec := List([1..RelRedLen],x->0);
                            for w in v do
                                RelReduceVec[Position(RelReduceLett,p+q+r+s+t+w)] := 1;
                            od;
                            Append(RelReduceMat,[RelReduceVec]);
                        od;
                    fi;
                od;
            od;
        od;
    od;
od;
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do
                for t in Concatenation([Gen0],GensLett) do
                    for u in Concatenation([Gen0],GensLett) do
                        if (p+q+r+s+t+u)*GenDeg1to4 = 6 then
                            for v in CupRel2Lett do
                                RelReduceVec := List([1..RelRedLen],x->0);
                                for w in v do
                                    RelReduceVec[Position(RelReduceLett,p+q+r+s+t+u+w)] := 1;
                                od;
                                Append(RelReduceMat,[RelReduceVec]);
                            od;
                        fi;
                    od;
                od;
            od;
        od;
    od;
od;


####
#### End preparation for relation reduction


rk:=RankMatrix(RelReduceMat*Z(2));



#Step-1 begins here: degree-7 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase7 do

    #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 1: ####
    ####
    uCocycle:=CB[7].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,7,1);
    ww:=[];
    for i in [1..(R!.dimension(8))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        #cupped :=Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]);      #then calculate u-cup-v
        
        #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(8))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[8].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 2 ended ####

        Lett1 := CupBase7Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        if cupped = List([1..Cohomology(TR,8)],x->0) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        else
            if CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase8,[cupped]);
                Append(CupBase8Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase8*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase8Lett) = false then
                    
            #### begin checking whether the relation is reducible from lower degree ones
            ####
            RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
            for w in solrel do
                #Print(w,"\n");
                solrelred := Position(RelReduceLett,w);
                            
                if solrelred = fail then      #if this relation contains new terms then we find a new relation
                    break;
                else
                    RelReduceVec[solrelred] := 1;
                fi;
            od;
            if (solrelred = fail) = false then
                rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                if (rk = rk1) = false then
                    solrelred :=fail;
                fi;
                #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
            fi;
            ####
            #### end checking whether the relation is reducible from lower degree ones

            if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new relation
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel8Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                Append(RelReduceMat,[RelReduceVec]);
                rk := rk1;
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-7 cup with degree-1-gen


#Step-2 begins here: (degree-3 cup with degree-3, or degree-2 cup with degree-2 cup with degree-2) cup with degree-2-gen
#
#
iu :=1;
for u in CupBase6 do
    if List([1..Length(Gen1)],x->CupBase6Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Ci-cup-Cj or Bi-cup-Bj-cup-Bk
        #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 1: ####
        ####
        uCocycle:=CB[6].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,6,2);
        ww:=[];
        for i in [1..(R!.dimension(8))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
            #cupped :=Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]);      #then calculate u-cup-v
        
            #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(8))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[8].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 2 ended ####
        

            Lett1 := CupBase6Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            if cupped = List([1..Cohomology(TR,8)],x->0) then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            else
                if CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else
                    sol :=SolutionMat(CupBase8*Z(2),cupped*Z(2));
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBase8,[cupped]);
                        Append(CupBase8Lett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase8Lett) = false then
               
                #### begin checking whether the relation is reducible from lower degree ones
                ####
                RelReduceVec := List([1..RelRedLen],x->0); #check whether this relation is reducible from lower ones
                for w in solrel do
                    #Print(w,"\n");
                    solrelred := Position(RelReduceLett,w);
                            
                    if solrelred = fail then      #if this relation contains new terms then we find a new relation
                        break;
                    else
                        RelReduceVec[solrelred] := 1;
                    fi;
                od;
                if (solrelred = fail) = false then
                    rk1:=RankMatrix(Concatenation(RelReduceMat*Z(2),[RelReduceVec]*Z(2)));
                    if (rk = rk1) = false then
                        solrelred :=fail;
                    fi;
                    #solrelred := SolutionMat(RelReduceMat*Z(2),RelReduceVec*Z(2));
                fi;
                ####
                #### end checking whether the relation is reducible from lower degree ones

                if solrelred = fail then          #if not reducible from lower degree ones, then we have found a new     relation
                    Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                    Append(CupRel8Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    Append(RelReduceMat,[RelReduceVec]);
                    rk := rk1;
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-3 cup with degree-3, degree-2 cup with degree-2 cup with degree-2) cup with degree-2-gen




#if CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
#    if (CupTemp = []) = false then
#        if (CupTemp[1] = []) = false then  #if cohomology dimension is not zero
#            Append(CupBase8,[CupTemp[1]]);
#            Append(CupBase8Lett,[CupTempLett[1]]);
#        fi;
#    fi;
#fi;


#if (CupBase8 = []) = false then

#    for cupped in CupTemp do                #check the cocycle-ness of those u-cup-v's that contain patterns in CupRelsLett
#        sol :=SolutionMat(CupBase8*Z(2),cupped*Z(2));
#        if sol = fail then                  #if u-cup-v is a genuine new cocycle
#            Append(CupBase8,[cupped]);
#            Append(CupBase8Lett,[Lett1]);
#        fi;
#    od;
#fi;



#Step-3 begins here: degree-4-gen cup with degree-4-gen
#
#
iu :=1;
for u in Gen4 do

    #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,4);
    ww:=[];
    for i in [1..(R!.dimension(8))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8])  -- part 1 ended ####

    iv :=1;
    for v in Gen4 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iu] + GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
            #IO := false;
            #for Lett2 in CupRelsLett do
            #    if NonNegativeVec(Lett1-Lett2) then #u-cup-v contains patterns in CupRelsLett
            #    IO :=true;       #then rewrite IO
            #    fi;
            #od;
        
            #cupped :=Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]);      #then calculate u-cup-v
            
            #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 2: ####
            ####
            vCocycle:=CB[4].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(8))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[8].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 2 ended ####
        
            if cupped = List([1..Cohomology(TR,8)],x->0) then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel8Lett,[[Lett1]]);
            elif CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase8,[cupped]);
                Append(CupBase8Lett,[Lett1]);
            else
                sol :=SolutionMat(CupBase8*Z(2),cupped*Z(2));
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                    if (Lett1 in CupBase8Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel8Lett,[solrel]);
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
#Step-3 finished: degree-4-gen cup with degree-4-gen

#Append(CupBase8,Gen8);
#Note that there's no Gen8!!
#for iu in [(1+Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7))..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7)+Length(Gen8))] do
#    if CupBase8 = [] then
#        Append(CupBase8,[Gen8[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7))]]);
#        Append(CupBase8Lett,[GensLett[iu]]);
#    else
#        sol :=SolutionMat(CupBase7*Z(2),Gen8[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7))]*Z(2));
#        if sol = fail then
#            Append(CupBase8,[Gen8[iu-(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6)+Length(Gen7))]]);
#            #Append(CupBase8Lett,[GensLett[iu]]);
#        else
#            Print("Error: containing fake degree-8 generator(s)!!\n");
#        fi;
#    fi;
#    #Append(CupBase8Lett,[GensLett[iu]]);
#od;



#Print("Chosen basis at degree 8:\n");
#PrintMonomialString(CupBase8Lett,GenDim1to4,",");

if Length(CupBase8) = Cohomology(TR,8) then
    Print("");#Print("dim(H^8)=", Cohomology(TR,8),".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^8) = ", Length(CupBase8) - Cohomology(TR,8),"\n");
fi;


fi;
#the above "fi;" is the one that ends working on the r=7 and r=8 relations, which is executed only when Gen4 is non-empty!
#
#
#
#
#####################################################################################

#
####   Begin printing cohomology ring   ####
#
mono:=List(List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4))],x->GensLett[x]),x->Letter2Monomial(x,GenDim1to4,GENNAMES[IT]));
Print("Mod-2 Cohomology Ring of Group No. ", IT, ":\n");
if (IT in [219,226]) then
    Print("Z2[", JoinStringsWithSeparator(mono,","), ",F1,F2]");
elif IT = 228 then
    Print("Z2[", JoinStringsWithSeparator(mono,","), ",F]");
else
    Print("Z2[", JoinStringsWithSeparator(mono,","), "]");
fi;



mono := 0;
if Length(CupRel2Lett)+Length(CupRel3Lett)+Length(CupRel4Lett)+Length(CupRel5Lett)+Length(CupRel6Lett)> 0 then
    Print("/<");
    if Length(CupRel2Lett)>0 then
        Print("R2");
        mono := 1;
    fi;
    if Length(CupRel3Lett)>0 then
        if mono = 1 then           #mono = 1 means at least one relation already printed
            Print(",R3");
        else
            Print("R3");
        fi;
        mono := 1;
    fi;
    if Length(CupRel4Lett)>0 then
        if mono = 1 then           #mono = 1 means at least one relation already printed
            Print(",R4");
        else
            Print("R4");
        fi;
        mono := 1;
    fi;
    if Length(CupRel5Lett)>0 then
        if mono = 1 then           #mono = 1 means at least one relation already printed
            Print(",R5");
        else
            Print("R5");
        fi;
        mono := 1;
    fi;
    if Length(CupRel6Lett)>0 then
        if mono = 1 then           #mono = 1 means at least one relation already printed
            Print(",R6");
        else
            Print("R6");
        fi;
    fi;
        Print(">\n");
fi;


if Length(CupRel2Lett) > 0 then
    #Print(Length(CupRel2Lett)," at deg 2: ");List(CupRel2Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R2:  ");List(CupRel2Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
if Length(CupRel3Lett) > 0 then
    #Print(Length(CupRel3Lett)," at deg 3: ");List(CupRel3Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R3:  ");List(CupRel3Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
if Length(CupRel4Lett) > 0 then
    #Print(Length(CupRel4Lett)," at deg 4: ");List(CupRel4Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R4:  ");List(CupRel4Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
if Length(CupRel5Lett) >0 then
    #Print(Length(CupRel5Lett)," at deg 5: ");List(CupRel5Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R5:  ");List(CupRel5Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;
if Length(CupRel6Lett) >0 then
    #Print(Length(CupRel6Lett)," at deg 6: ");List(CupRel6Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    Print("R6:  ");List(CupRel6Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
fi;

#Begin printing Degree 7 and 8 relations:
#
#
if Length(Gen4) >0 and Length(Size(R)) >= 9 then
    if Length(CupRel7Lett) >0 then
        #Print(Length(CupRel7Lett)," at deg 7: ");List(CupRel7Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
        Print("R7:  ");List(CupRel7Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    fi;
    if Length(CupRel8Lett) >0 then
        #Print(Length(CupRel8Lett)," at deg 8: ");List(CupRel8Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
        Print("R8:  ");List(CupRel8Lett,x->PrintMonomialString(x,GenDim1to4,"+",GENNAMES[IT]));Print("\n");
    fi;
fi;
#
#
#End printing Degree 7 and 8 relations:

####   End printing cohomology ring   ####
#


#Print("End Group No. ", IT, ".\n");
Print("===========================================\n");


return [CupBase1Lett,CupBase2Lett,CupBase3Lett,CupBase4Lett];
#rec(GensAtDegN:=List([1..4],x->List([1..Length(Gens[x])],y->Letters[x,y])),RelsAtDegN:=[CupRel2Lett,CupRel3Lett,CupRel4Lett,CupRel5Lett,CupRel6Lett]);
end;
#####################################################################
#####################################################################



#####################################################################
#####################################################################

PointGroupTranslationExtension:=function(arg)
local
    Gs,arithmeticNo,ZZZ,Gpt,R,C,Homo,
    T1,T2,T3,elem,eleml,conjT1,conjT2,conjT3,
    i,j,k,l,rep,rep1,flag;
    
#Z1:=Group([()]);
#Z2:=Group([(1,2)]);
#Z2Z2:=Group([(1,2),(3,4)]);
#Z2Z2Z2:=Group([(1,2),(3,4),(5,6)]);
#Z4:=Group([(1,2),(1,2,3,4)]);
#Z4Z2:=Group([(1,2),(1,2,3,4),(5,6)]);
#Dih4:=Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6)]);
#Dih4Z2:=Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6),(7,8)]);
#Z3:=Group([(1,2,3)]);
#Z3Z2:=Group([(1,2,3),(4,5)]);
#Dih3:=Group([(1,2,3),(2,3)(4,5)]);
#Dih3Z2:=Group([(1,2,3),(2,3)(4,5),(6,7)]);
#Z3Z2Z2:=Group([(1,2,3),(4,5),(6,7)]);
#Dih3Z2Z2:=Group([(1,2,3),(2,3)(4,5),(6,7),(8,9)]);
#A4:=Group([(1,2)(3,4),(1,3)(2,4),(1,2,3)]);
#A4Z2:=Group([(1,2)(3,4),(1,3)(2,4),(1,2,3),(5,6)]);
#S4:=Group([(1,2)(3,4),(1,3)(2,4),(1,2,3),(1,2)]);
#S4Z2:=Group([(1,2)(3,4),(1,3)(2,4),(1,2,3),(1,2),(5,6)]);

Gs:=[Group([(1,2)]),Group([(1,2),(3,4)]),Group([(1,2),(3,4),(5,6)]),Group([(1,3)(2,4),(1,2,3,4)]),Group([(1,3)(2,4),(1,2,3,4),(5,6)]),Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6)]),Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6),(7,8)]),Group([(1,2,3)]),Group([(1,2,3),(4,5)]),Group([(1,2,3),(2,3)(4,5)]),Group([(1,2,3),(2,3)(4,5),(6,7)]),Group([(1,2,3),(4,5),(6,7)]),Group([(1,2,3),(2,3)(4,5),(6,7),(8,9)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(5,6)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(1,2)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(1,2),(5,6)])];

arithmeticNo:=[[[2],[3,4],[5],[6,7],[8,9]],[[10,11,13,14],[12,15],[16,17,18,19],[20,21],[22],[23,24],[25,26,27,28,29,30,31,32,33,34],[35,36,37],[38,39,40,41],[42,43],[44,45,46]],[[47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62],[63,64,65,66,67,68],[69,70],[71,72,73,74]],[[75,76,77,78],[79,80],[81],[82]],[[83,84,85,86],[87,88]],[[89,90,91,92,93,94,95,96],[97,98],[99,100,101,102,103,104,105,106],[107,108,109,110],[111,112,113,114],[115,116,117,118],[119,120],[121,122]],[[123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138],[139,140,141,142]],[[143,144,145],[146]],[[147],[148],[168,169,170,171,172,173],[174]],[[149],[150],[151],[152],[153],[154],[155],[156],[157],[158],[159],[160,161]],[[162,163],[164,165],[166,167],[177,178,179,180,181,182],[183,184,185,186],[187,188],[189,190]],[[175,176]],[[191,192,193,194]],[[195],[196],[197],[198],[199]],[[200,201],[202,203],[204],[205],[206]],[[207,208],[209,210],[211],[212,213],[214],[215],[216],[217],[218],[219],[220]],[[221,222,223,224],[225,226,227,228],[229,230]]];

ZZZ:=GL(3,Integers);;

T1:=[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]]; #standard translation T1
T2:=[[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]]; #standard translation T2
T3:=[[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]]; #standard translation T3

Print("Group Extension Information for the 73 Arithmetic Classes Carrying Representation rho:\n");
Print("H^2_rho(PG,Z^3):\n");


for i in [1..Length(arithmeticNo)] do           #Length(arithmeticNo) == 17
    Gpt:=Gs[i];
    R:=ResolutionFiniteGroup(Gpt,4);
    for j in [1..Length(arithmeticNo[i])] do
        flag := arithmeticNo[i][j][1];
        
        for k in arithmeticNo[i][j] do
            for l in [1..Length(PGGens230[k])] do
                elem := List([1..3],x->List([1..3],y->PGGens230[flag][l][x][y]));
                eleml := List([1..3],x->List([1..3],y->PGGens230[k][l][x][y]));
                if (elem = eleml) = false then
                    Print("point group generators not chosen consistently with others in the same arithmetic class!\n");
                fi;
            od;
        od;
        
        rep := [];
        for elem in PGGens230[flag] do
            conjT1 := elem * T1 * elem^(-1);
            conjT2 := elem * T2 * elem^(-1);
            conjT3 := elem * T3 * elem^(-1);
            Append(rep,[[[conjT1[1,4],conjT2[1,4],conjT3[1,4]],[conjT1[2,4],conjT2[2,4],conjT3[2,4]],[conjT1[3,4],conjT2[3,4],conjT3[3,4]]]]);
        od;
        
        Homo := GroupHomomorphismByImages(Gpt,ZZZ,GeneratorsOfGroup(Gpt),rep);
        C:=HomToIntegralModule(R,Homo);
        if Length(arithmeticNo[i][j]) = 1 then
            Print("No. ", arithmeticNo[i][j][1],": ", Cohomology(C,2),"  ");
        elif List([1..Length(arithmeticNo[i][j])],x->arithmeticNo[i][j][x]) = List([1..Length(arithmeticNo[i][j])],x->arithmeticNo[i][j][1]+x-1) then
            Print("No. ", arithmeticNo[i][j][1],"-", arithmeticNo[i][j][Length(arithmeticNo[i][j])],": ", Cohomology(C,2),"  ");
        else
            Print("No. ", JoinStringsWithSeparator(List([1..Length(arithmeticNo[i][j])],x->String(arithmeticNo[i][j][x])),"&"),": ", Cohomology(C,2),"  ");
        fi;
    od;
    Print("\n");
od;


return true;
end;
#####################################################################
#####################################################################




#####################################################################
#####################################################################
IrreducibleWyckoffPoints:=function(arg)
local
    IT, Dim, SG, Rec, IWPs, x,
    WyckoffPosRelations, WyckoffGraphRecord;
    
#####################################################################
WyckoffPosRelations := function( W )
    local S, T, d, len, gens, G, L, m, O, o, i, j, k, Si, Sj, index, lst;

    S := WyckoffSpaceGroup( W[1] );
    T := TranslationBasis( S );
    d := DimensionOfMatrixGroup( S ) - 1;
    len  := Length( W );
    gens := GeneratorsOfGroup( S );
    gens := Filtered( gens, g -> g{[1..d]}{[1..d]} <> IdentityMat( d ) );
    if IsAffineCrystGroupOnLeft( S ) then
        gens := List( gens, TransposedMat );
    fi;
    G := GroupByGenerators( gens, One( S ) );
    L := List( W, w -> rec( translation := WyckoffTranslation( w ),
                            basis       := WyckoffBasis( w ),
                            spaceGroup  := S ) );

    m := NullMat( len, len );
    for i in [1..len] do
        O := Orbit( G, L[i], ImageAffineSubspaceLattice );
        for j in [1..len] do
            Sj := WyckoffStabilizer( W[j] );
            Si := WyckoffStabilizer( W[i] );
            index := Size(Sj) / Size(Si);
            if Length(L[j].basis) < Length(L[i].basis) and IsInt(index) then
                lst := Filtered(O,o->IsSubspaceAffineSubspaceLattice(o,L[j]));
                m[j][i] := Length( lst );
            fi;
        od;
    od;

    for i in Reversed([1..Length(W)]) do
        for j in Reversed([1..i-1]) do
            if m[j][i]<>0 then
                for k in [1..j-1] do
                    if m[k][j]<>0 then m[k][i]:=0; fi;
                od;
            fi;
        od;
    od;

    return m;

end;

#############################################################################
##
#F  WyckoffGraphRecord( <lst> ) . . . . . . . Create record for Wyckoff graph
##
WyckoffGraphRecord := function( lst )

    local L, m, R, i, level, j;

    L := List( lst, w -> rec( wypos := w,
                              dim   := Length( WyckoffBasis(w) ),
                              size  := Size( WyckoffStabilizer(w) ),
                              class := w!.class ) );
    Sort( L, function(a,b) return a.size > b.size; end );

    m := WyckoffPosRelations( List( L, x -> x.wypos ) );

    R := rec( levels   := [],
              classes  := [],
              vertices := [],
              edges    := [] );

    for i in [1..Length(L)] do
        level := [ L[i].dim, L[i].size ];
        AddSet( R.levels, level );
        AddSet( R.classes, [ L[i].class, level ] );
        Add( R.vertices, [ L[i].wypos, level, L[i].class ] );
        for j in [1..i-1] do
            if m[j][i]<>0 then Add( R.edges, [ i, j, m[j][i] ] ); fi;
        od;
    od;

    return R;

end;
#####################################################################

IT := arg[1];
if Length(arg) = 1 then
    SG := SpaceGroupIT(3,IT);
else
    SG := SpaceGroupIT(arg[2],IT);  #arg[2] must give the dimension (either 2 or 3);
fi;

Rec := WyckoffGraphRecord(WyckoffPositions(SG));


IWPs := Difference([1..Length(Rec.vertices)],Set(List(Rec.edges,x->x[1]))); #IWPs stores the index of the IWPs;

for x in IWPs do
    Print(Rec.vertices[x][1],": ", ["point","line","plane","volume"][1 + Rec.vertices[x][2][1]],"\n");
od;

Print("Number of IWPs for Group No.", IT, ": ",Length(IWPs),"\n");

return Length(IWPs);
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
    GenName_standard, Rdim,lst1,lst2,lst3,lst4,lst5,lst6,lst1to4,
    MatToPow,GapToPow,Fbarhomotopyindinv,Invofg,Prodg1g2Pow,IndToElem,
    Homotopydeg1,Homotopydeg2,Homotopydeg3,Homotopydeg4,
    func,funcs,receive,FuncVal,TopoInvdeg3,
    Gen1, Gen2, Gen3, Gen4, GensGAP, GensDim1to4, GensDeg1to4,
    BasesLett, Base1Lett, Base2Lett, Base3Lett, Base4Lett, LSMLett,
    g1,g2,g3,overcomplete_g,Mat,mat1,mat2,vec,sol,LSMMat,CountLSM,
    v1,v2,v3,x1,y1,z1,x2,y2,z2,x3,y3,z3,
    i,j,k,p,x,y;

#####################################################################
IndToElem:=function(lst,Lst)         #given lst of 0&1's, extract the elements in Lst corresponding to 1's
local i,l,out;

if (Length(lst) = Length(Lst)) = false then
    Print("IndToElem:  Input is wrong!!\n");
else
    l:=Length(lst);
    out := [];
    for i in [1..l] do
        if lst[i] = 1 then
            Append(out,[Lst[i]]);
        fi;
    od;
fi;
return out;
end;
#####################################################################
MatToPow:=function(mat)            #given 4x4 matrix, output the list of powers of group generators
local i, mat33, trans;

mat33:=List([1..3],i->List([1..3],j->mat[i,j]));

i:=Position(PGMat33,mat33);

trans:=mat*PGMatinv[i];

return Concatenation(List([1..3],x->trans[x,4]),PGind[i]);
end;
#####################################################################
GapToPow:=function(i)            #given the index i s.t. mat:=R!.elts[i], output the list of powers of group generators

return MatToPow(TransposedMat(PreImage(Gp,R!.elts[i])));
end;
#####################################################################
Invofg:=function(v)
local vpg, transmat;

transmat := [[1,0,0,v[1]],[0,1,0,v[2]],[0,0,1,v[3]],[0,0,0,1]];
vpg := List([4..(Length(v))],x->v[x]);

return MatToPow((transmat * PGMatinv[Position(PGind,vpg)]^(-1))^(-1));
end;
#####################################################################
Prodg1g2Pow:=function(v1,v2)
local vpg1, vpg2, transmat1, transmat2, prod;

transmat1 := [[1,0,0,v1[1]],[0,1,0,v1[2]],[0,0,1,v1[3]],[0,0,0,1]];
transmat2 := [[1,0,0,v2[1]],[0,1,0,v2[2]],[0,0,1,v2[3]],[0,0,0,1]];
vpg1 := List([4..(Length(v1))],x->v1[x]);
vpg2 := List([4..(Length(v2))],x->v2[x]);

prod := transmat1 * PGMatinv[Position(PGind,vpg1)]^(-1) * transmat2 * PGMatinv[Position(PGind,vpg2)]^(-1);

return MatToPow(prod);
end;
#####################################################################
Fbarhomotopyindinv:=function(i,lst)            #This is the function Fbarhomotopyindinv in Mathematica

return List(lst,x->Concatenation(x,[i]));
end;
#####################################################################
FuncVal:=function(lett,v)                #Given a monomial of degree 1, 2, or 3, and argument (for the degree 3 monomial, the argument is either g1,g1,g1 or g1,g2,g2 or g1,g2,g3), evaluate the cocycle.
local ct,deg,i,ival,j,jval,k,val,lett1;
deg := GensDeg1to4*lett;

if deg = 0 then
    Print("Degree of gen is wrong!!\n");
fi;


lett1 := ShallowCopy(lett);

for i in [1..Length(lett)] do            #finding the first generator that exists in lett
    if lett1[i] > 0 then
        ival := lett1[i];                #label of generator stored in i; power of this generator stored in ival
        break;
    fi;
od;

if deg = 1 then                          #if evaluating a 1-cocycle
    val := funcs[1][i](v[1]);
elif deg = 2 then                        #if evaluating a 2-cocycle
    if i > GensDim1to4[1] then
    val := funcs[2][i-GensDim1to4[1]](v[1],v[2]);
    else
        lett1[i] := lett1[i] - 1;
        j := Position(lett1,1);
        val := funcs[1][i](v[1]) * funcs[2][j](v[2]);
    fi;
elif deg = 3 then                        #if evaluating a 3-cocycle
    if i > GensDim1to4[1]+GensDim1to4[2] then                              #if a degree-3 generator
        val := funcs[3][i-GensDim1to4[1]-GensDim1to4[2]](v[1],v[2],v[3]);
    else                                                                   #if not a degree-3 gen., then must be a cup prod.
        lett1[i] := lett1[i] - 1;
        for j in [1..Length(lett)] do
            if lett1[j] > 0 then
                jval := lett1[j];        #label of generator stored in j; power of this generator stored in jval
                break;
            fi;
        od;
        if j > GensDim1to4[1] then                                         #if a degree-1 gen. cup a degree-2 gen.
            val := funcs[1][i](v[1]) * funcs[2][j-GensDim1to4[1]](v[2],v[3]);
        else                                                               #if not, then must be cup of three degree-1 gens.
            lett1[j] := lett1[j] - 1;
            k := Position(lett1,1);
            val := funcs[1][i](v[1]) * funcs[1][j](v[2]) * funcs[1][k](v[3]);
        fi;
    fi;
else
    Print("Error: Evaluating 4- or higher cocycles not implemented yet!! \n");
fi;
return val;
end;
#####################################################################
TopoInvdeg3:=function(arg) #usage: TopoInvdeg3(list_of_group_elements,list_of_letters,[matrices giving linear combination of letters])
local gs,letters,solrels, vallist;

gs := arg[1];               #List of group elements [g1] or [g1,g2] or [g1,g2,g3] at which the cocycles are evaluated
letters := arg[2];          #List of letters representing the monomials


if Length(arg) = 2 then
    solrels := List([1..Length(letters)],x->List([1..Length(letters)],y->0));
    for i in [1..Length(letters)] do
        solrels[i][i] := 1;
    od;
else
    solrels := arg[3];
fi;

if Length(gs) = 1 then
    if Prodg1g2Pow(gs[1],gs[1]) = gs[1]*0 then
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[1],gs[1]]));                                    #Topo inv varphi1
    else
        Print("varphi(g,g,g) is not a topological invariant!!!!\n");
    fi;
elif Length(gs) = 2 then
    if (Prodg1g2Pow(gs[2],gs[2]) = gs[2]*0) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) then #Topo inv varphi2
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],gs[2]])+FuncVal(x,[gs[2],gs[1],gs[2]])+FuncVal(x,[gs[2],gs[2],gs[1]]));
    else
        Print("varphi(g1,g2,g2) is not a topological invariant!!!!\n");
        Print(gs[1],gs[2],"\n");
    fi;

elif Length(gs) = 3 then

    if (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[1],gs[3]) = Prodg1g2Pow(gs[3],gs[1])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2])) then     #Topo inv varphi3
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],gs[3]])+FuncVal(x,[gs[1],gs[3],gs[2]])+FuncVal(x,[gs[2],gs[1],gs[3]])+FuncVal(x,[gs[2],gs[3],gs[1]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));
    
    elif ((Prodg1g2Pow(gs[2],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[2])) and (Prodg1g2Pow(gs[1],gs[3]) = Prodg1g2Pow(gs[3],gs[1])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2]))) then     #Topo inv tildevarphi for No. 7,26,36,39,46,57,62
        
        vallist := List(letters,x->FuncVal(x,[gs[3],Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[1],Invofg(gs[2]))])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[1],Invofg(gs[2])])+FuncVal(x,[gs[3],gs[2],Invofg(gs[2])])+FuncVal(x,[Prodg1g2Pow(gs[1],gs[2]),gs[3],Prodg1g2Pow(gs[1],Invofg(gs[2]))])+FuncVal(x,[gs[1],gs[3],gs[2]])+FuncVal(x,[gs[1],gs[3],Invofg(gs[2])])+FuncVal(x,[gs[2],gs[3],Invofg(gs[2])])+FuncVal(x,[Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[1],Invofg(gs[2])),gs[3]])+FuncVal(x,[gs[1],gs[2],gs[3]])+FuncVal(x,[gs[1],Invofg(gs[2]),gs[3]])+FuncVal(x,[gs[2],Invofg(gs[2]),gs[3]]));

    elif ((Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(gs[2],gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 76 & 78
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Invofg(gs[1]),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[1]])+FuncVal(x,[gs[2],gs[3],gs[2]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));

    elif ((Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 4
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));

    elif ((Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 9,161
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[1]])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[2]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));
    else
        Print("varphi(g1,g2,g3) is not a topological invariant!!!!\n");
    fi;

elif Length(gs) = 4 then
    if (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[3],gs[3]) = gs[2]) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) then     #Topo inv varphi3 for No. 19 and 198
        
        vallist := List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])]) + FuncVal(x,[gs[1],gs[2],Invofg(gs[1])]) + FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]]) + FuncVal(x,[gs[1],gs[2],Invofg(gs[2])])+ FuncVal(x,[gs[2],gs[1],Invofg(gs[2])]) + FuncVal(x,[gs[2],Invofg(gs[2]),gs[1]]) + FuncVal(x,[gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1]]) + FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]]));
    
    elif (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) and (Prodg1g2Pow(gs[4],gs[3]) = Prodg1g2Pow(Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[3],gs[4]))) then     #Topo inv varphi3 for No. 29

        vallist := List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])])+FuncVal(x,[gs[1],gs[2],Invofg(gs[1])])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[3],Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[3],gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[3],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]));
    
    elif (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) and (Prodg1g2Pow(gs[4],gs[3]) = Prodg1g2Pow(gs[1],Prodg1g2Pow(gs[3],gs[4]))) then     #Topo inv varphi3 for No. 33
        
        vallist :=
        List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])])+FuncVal(x,[gs[1],gs[2],Invofg(gs[1])])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[2],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[1],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[3],Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[3],gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[3],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]));
    fi;
else
    Print("Wrong in checking topological invariant: Number of group elements is not between 1 and 4!!\n");
fi;
return GF2ToZ((solrels*vallist)*Z(2));
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
#Read("~/Downloads/Space_Group_Cocycles.gi");


if 1<=IT and IT<=230 then
    PGGen:=PGGens230[IT];
    funcs:=funcs230[IT];
else
    Print("Space Group IT not found!!!!", "\n");
fi;

PGGen33 := List([1..Length(PGGen)],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));

if Length(PGGen) = 0 then
    Append(PGind,[[]]);
    PGMat33:=[[[1,0,0],[0,1,0],[0,0,1]]];
    PGMatinv:=[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]];
elif Length(PGGen) = 1 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        Append(PGind,[[o1]]);
        Append(PGMat33,[PGGen33[1]^o1]);
        Append(PGMatinv,[(PGGen[1]^o1)^(-1)]);
    od;
elif Length(PGGen) = 2 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #the ordering of the group generators must strictly follow this condition
            Append(PGind,[[o1,o2]]);
            Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2]);
            Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2)^(-1)]);
        od;
    od;
elif Length(PGGen) = 3 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #the ordering of the group generators must strictly follow this condition
            for o3 in [0..(Order(PGGen33[3])-1)] do
                Append(PGind,[[o1,o2,o3]]);
                Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3]);
                Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3)^(-1)]);
            od;
        od;
    od;
elif Length(PGGen) = 4 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #the ordering of the group generators must strictly follow this condition
            for o3 in [0..(Order(PGGen33[3])-1)] do
                for o4 in [0..(Order(PGGen33[4])-1)] do
                    Append(PGind,[[o1,o2,o3,o4]]);
                    Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3*PGGen33[4]^o4]);
                    Append(PGMatinv,[(PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3*PGGen[4]^o4)^(-1)]);
                od;
            od;
        od;
    od;
elif Length(PGGen) = 5 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #the ordering of the group generators must strictly follow this condition
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
else
    Print("Number of Point Group Generators Exceeds 5 -- WRONG!!!");
fi;



#Standard resolution is generated by the following:
G:= Group(Concatenation([T1,T2,T3],PGGen));
Gp:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G))));

#Non-standard resolution is generated by the following:
#Gp:=IsomorphismPcpGroup(SpaceGroupBBNWZ(3,IT));
if IT=222 then
    Gp:=IsomorphismPcpGroup(SpaceGroupIT(3,222));
fi;

#Below constructs the resolution for the group:
#
if (IT in [108, 109, 120, 130, 136, 140, 142, 197, 204]) = true then
    R:=ResolutionAlmostCrystalGroup(Image(Gp),9);
elif IT <= 220 then
    R:=ResolutionAlmostCrystalGroup(Image(Gp),7);
else
    R:=ResolutionAlmostCrystalGroup(Image(Gp),6);
fi;
#
#Resolution for the group now constructed.



Homotopydeg1:=List([1..R!.dimension(1)],x->List(R!.boundary(1,x),y->[y[2]]));
Homotopydeg2:=List([1..R!.dimension(2)],x->Concatenation(List(R!.boundary(2,x),y->Fbarhomotopyindinv(y[2],Homotopydeg1[AbsInt(y[1])]))));
Homotopydeg3:=List([1..R!.dimension(3)],x->Concatenation(List(R!.boundary(3,x),y->Fbarhomotopyindinv(y[2],Homotopydeg2[AbsInt(y[1])]))));


CB:=[];
for p in [1..4] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

#Print(IT,":  \n");
#Print(List([1..4],x->Cohomology(HomToIntegersModP(R,2),x),"\n");
#Print(CB[1].cocyclesBasis);
#Print("\n");



Gen1:=[];
for func in funcs[1] do
    Append(Gen1,[CB[1].cocycleToClass(List([1..R!.dimension(1)],x->RemInt(Sum(List(Homotopydeg1[x],y->func(GapToPow(y[1])))),2)))]);
od;


Gen2:=[];
for func in funcs[2] do
    Append(Gen2,[CB[2].cocycleToClass(List([1..R!.dimension(2)],x->RemInt(Sum(List(Homotopydeg2[x],y->func(GapToPow(y[1]),GapToPow(y[2])))),2)))]);
od;

Gen3:=[];

if (IT in [225,227,229]) = true then
    Gen3[1] := CB[3].cocycleToClass(List([1..R!.dimension(3)],x->RemInt(Sum(List(Homotopydeg3[x],y->funcs[3][1](GapToPow(y[1]),GapToPow(y[2]),GapToPow(y[3])))),2)));
    if IT = 225 then
        Gen3[2] := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];    #Must use standard resolution!!
        Gen3[3] := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];    #Must use standard resolution!!
    elif IT = 227 then
        Gen3[2] := [0, 0, 0, 0, 0, 0, 1, 0, 0];                #Must use standard resolution!!
    elif IT = 229 then
        Gen3[2] := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];    #Must use standard resolution!!
    fi;
else
    for func in funcs[3] do
        Append(Gen3,[CB[3].cocycleToClass(List([1..R!.dimension(3)],x->RemInt(Sum(List(Homotopydeg3[x],y->func(GapToPow(y[1]),GapToPow(y[2]),GapToPow(y[3])))),2)))]);
    od;
fi;


Gen4:=[];


GensGAP:=Mod2RingGenerators(R,4,3); #GensGAP: Generators of the mod-2 cohomology ring at degree 1-4 that are worked out by GAP; GensGAP is to be compared with those worked out from explicit cochain expressions (Gen1,Gen2,Gen3). Here: R is the resolution, 4 is the maximal degree at which the generators are worked out, and 3 is the dimension of the space group (if 2 instead then computes wallpaper group).

#Print("Matching generators for space group No. ", IT, "\n");
#Print("basis given by gap is:", GensGAP, "\n");
#Print("basis given by func is:", [Gen1,Gen2,Gen3], "\n");

#if (Length(Gen1) = Length(GensGAP[1])) = false then
#    Print("Number of Degree-1 generators does not match!!!:", Length(Gen1),"!=",Length(GensGAP[1]),"\n");
#elif (Length(Gen2) = Length(GensGAP[2])) = false then
#    Print("Number of Degree-2 generators does not match!!!:", Length(Gen2),"!=",Length(GensGAP[2]),"\n");
#elif (Length(Gen3) = Length(GensGAP[3])) = false then
#    Print("Number of Degree-3 generators does not match!!!:", Length(Gen3),"!=",Length(GensGAP[3]),"\n");
#else
#    #Print("Generators at degree 1,2,3 matched.","\n");
#    Print("\n");
#fi;


if (IT in [108, 109, 120, 130, 136, 140, 142, 197, 204, 230]) = true then
    Gen4 := GensGAP[4];
    if IT = 108 then
        Gen4 := [[0, 0, 1, 1, 1, 0, 1]];                                                  #Must use standard resolution!!
    elif IT = 120 then
        Gen4 := [[0, 1, 1, 0, 0, 1, 1, 0],[0, 0, 0, 0, 0, 1, 0, 1]];                      #Must use standard resolution!!
    elif IT = 140 then
        Gen4 := [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1]];     #Must use standard resolution!!
    elif IT = 142 then
        Gen4 := [[0, 1, 1, 0, 0, 0, 1]];                                                  #Must use standard resolution!!
    elif IT = 230 then
        Gen4 := [[0, 1, 1, 1, 1]];                                                        #Must use standard resolution!!
    fi;
fi;

#Gen1:=GensGAP[1];
#Gen2:=GensGAP[2];
#Gen3:=GensGAP[3];


################## BELOW ARE CODES RELATED TO NAMING OF COCYCLE GENERATORS ######################

#
#### Suggested name replacement:####
#

#Find the highest 1 in the cocycle vector to infer the place in the standard translation-point group LHS spectral sequence,
#And use this to label the cocycles.
#2-cocycles at (0,2),(1,1),(2,0) are labeled as B_{translation-type}, B_beta, B_alpha;
#3-cocycles at (0,3),(1,2),(2,1),(3,0) are labeled by C_{xyz}, C_gamma, C_beta, C_alpha;
#4-cocycles at (1,3),(2,2),(3,1),(4,0) are labeled by D_delta, D_gamma, D_beta, D_alpha.
#and so on.


if false then


Rdim := List([1..4],x->R!.dimension(x));
lst1 := [Rdim[1]-3];
lst2 := [Rdim[2]-3-3*lst1[1],Rdim[2]-3];
lst3 := [Rdim[3]-1-3*(lst1[1]+lst2[1]),Rdim[3]-1-3*lst1[1],Rdim[3]-1];
lst4 := [Rdim[4]-lst1[1]-3*(lst2[1]+lst3[1]),Rdim[4]-lst1[1]-3*lst2[1],Rdim[4]-lst1[1]];
#lst5 := [Rdim[5]-lst2[1]-3*(lst3[1]+lst4[1]),Rdim[5]-lst2[1]-3*lst3[1],Rdim[5]-lst2[1]];
#lst6 := [Rdim[6]-lst3[1]-3*(lst4[1]+lst5[1]),Rdim[6]-lst3[1]-3*lst4[1],Rdim[6]-lst3[1]];
lst1to4 := [lst1,lst2,lst3,lst4];

GenName_standard := [];
#GenName_standard := List([1..Length(Gen1)],x->JoinStringsWithSeparator([GENNAMES[IT][x],"in",IT],""));

j := Length(Gen1)+1;
for i in [2..4] do
    for o1 in [Gen1,Gen2,Gen3,Gen4][i] do
        v1 := Positions(CB[i].classToCocycle(o1),1);
        v2 := v1[Length(v1)];
        if v2 <=lst1to4[i][1] then
            Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"alphain",IT," :=",GENNAMES[IT][j],"in",IT],"")]);
            #Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"alphain",IT],"")]);
        elif v2<=lst1to4[i][2] then
            Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"betain",IT," :=",GENNAMES[IT][j],"in",IT],"")]);
            #Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"betain",IT],"")]);
        elif i=2 then
            o2:=CB[i].classToCocycle(o1);
            o3 := [o2[Length(o2)-2],o2[Length(o2)-1],o2[Length(o2)]];
            if o3 = [1,0,0] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Bxyin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Bxyin", IT],"")]);
            elif o3 = [0,0,1] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Byzin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Byzin", IT],"")]);
            elif o3 = [0,1,1] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Bzxyin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Bzxyin", IT],"")]);
            elif o3 = [1,0,1] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Byxzin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Byxzin", IT],"")]);
            elif o3 = [1,1,0] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Bxyzin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Bxyzin", IT],"")]);
            elif o3 = [1,1,1] then
                Append(GenName_standard,[JoinStringsWithSeparator(["Bxyxzyzin", IT, " :=",GENNAMES[IT][j],"in",IT],"")]);
                #Append(GenName_standard,[JoinStringsWithSeparator(["Bxyxzyzin", IT],"")]);
            fi;
        elif v2<=lst1to4[i][3] then
            Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"gammain",IT," :=",GENNAMES[IT][j],"in",IT],"")]);
            #Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"gammain",IT],"")]);
        elif i=3 then
            Append(GenName_standard,[JoinStringsWithSeparator(["Cxyzin",IT," :=",GENNAMES[IT][j],"in",IT],"")]);
            #Append(GenName_standard,[JoinStringsWithSeparator(["Cxyzin",IT],"")]);
        else # v2>lst1to4[i][3] and i>3
            Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"deltain",IT," :=",GENNAMES[IT][j],"in",IT],"")]);
            #Append(GenName_standard,[JoinStringsWithSeparator([["A","B","C","D"][i],"deltain",IT],"")]);
        fi;
        j :=j+1;
    od;
od;

Print("(*",IT,":*)",GenName_standard);

fi;

#
####
#
################## ABOVE ARE CODES RELATED TO NAMING OF COCYCLE GENERATORS ######################






BasesLett := Mod2RingGensAndRels(IT,3,R,[Gen1,Gen2,Gen3,Gen4]);

GensDim1to4 := [Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4)];
GensDeg1to4:=Concatenation(List([1..Length(Gen1)],x->1),List([1..Length(Gen2)],x->2),List([1..Length(Gen3)],x->3),List([1..Length(Gen4)],x->4));

Base1Lett := BasesLett[1];
Base2Lett := BasesLett[2];
Base3Lett := BasesLett[3];
Base4Lett := BasesLett[4];


#Print("Degree-4 generator:", GensGAP[4],"\n");
#Print("Degree-5 generator:", GensGAP[5],"\n");
#Print("Degree-6 generator:", GensGAP[6],"\n");


#Print all the elements of the mod-2 cohomology at degree 3:
#Print("List of degree-3 elements: \n");
#PrintMonomialString(Base3Lett,GensDim1to4,",",GENNAMES[IT]);




############################### BELOW ARE LSM RELATED CODES ###############################



overcomplete_g:=[];
Mat:=[];



#First: record all the LSM TIs, which have been given in Space_Group_Cocycles.gi
#
#
CountLSM := [];
for x in IWP[IT] do
    if (x[2] = []) = false then
        Append(Mat,[TopoInvdeg3(x[2],Base3Lett)]);
        Append(overcomplete_g,[x[2]]);
        Append(CountLSM,[x[2]]);
    fi;
od;

#Print("LSM topo invariants just added. Now the rank is: ", RankMatrix(Mat*Z(2)),"\n");
#Print(List(LSMMat.vectors,x->GF2ToZ(x)),"\n");



#Second: find all the non-LSM TIs, which are of one of the following four types:
#
#
for v2 in PGind do
    if (IT <= 220) or (v2[3] = 0) then
        for x2 in [-2..2] do
            for y2 in [-2..2] do
                for z2 in [-2..2] do
                    g2 := Concatenation([x2,y2,z2],v2);
                    mat2 := [[1,0,0,x2],[0,1,0,y2],[0,0,1,z2],[0,0,0,1]] * PGMatinv[Position(PGind,v2)]^(-1);
                    if (mat2^2 = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) then
                        if Trace(mat2)=0 then                    #C2 rotation
                            vec := TopoInvdeg3([g2],Base3Lett);
                            sol :=SolutionMat(Mat*Z(2),vec*Z(2));
                            if sol = fail then             #then we find a new non-LSM topo inv associated with C2 rotation
                                Append(Mat,[vec]);
                                Append(overcomplete_g,[[g2]]);
                            fi;
                        elif Trace(mat2)=2 then                  #Mirror
                            vec := TopoInvdeg3([g2],Base3Lett);
                            sol :=SolutionMat(Mat*Z(2),vec*Z(2));
                            if sol = fail then             #then we find a new non-LSM topo inv associated with mirror
                                Append(Mat,[vec]);
                                Append(overcomplete_g,[[g2]]);
                            fi;
                            for v1 in PGind do
                                for x1 in [-2..2] do
                                    for y1 in [-2..2] do
                                        for z1 in [-2..2] do
                                            g1 := Concatenation([x1,y1,z1],v1);
                                            if ((g1 = (g1*0)) = false) and Prodg1g2Pow(g2,g1) = Prodg1g2Pow(g1,g2) then
                                                vec := TopoInvdeg3([g1,g2],Base3Lett);
                                                sol :=SolutionMat(Mat*Z(2),vec*Z(2));
                                                if sol = fail then       #then we find a new non-LSM topo inv associated with   commuting couples g1 and g2, where g2 is a mirror
                                                    Append(Mat,[vec]);
                                                    Append(overcomplete_g,[[g1,g2]]);
                                                fi;
                                            fi;
                                        od;
                                    od;
                                od;
                            od;
                        fi;
                    elif (mat2^4 = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) and (Trace(mat2)=2) then         #C4 rotation
                        vec := TopoInvdeg3([g2,Prodg1g2Pow(g2,g2)],Base3Lett);
                        sol :=SolutionMat(Mat*Z(2),vec*Z(2));
                        if sol = fail then       #then we find a new non-LSM topo inv associated with commuting couples C4 and   C4^2
                            Append(Mat,[vec]);
                            Append(overcomplete_g,[[g2,Prodg1g2Pow(g2,g2)]]);
                        fi;
                    fi;
                od;
            od;
        od;
    fi;
od;


#Print(List(Mat*Z(2),x->GF2ToZ(x)));




if RankMatrix(Mat*Z(2)) = Length(Base3Lett) and RankMatrix(Mat*Z(2)) = Length(Mat) then

    #Print("Full Rank achieved: ", RankMatrix(Mat*Z(2)),"=", Length(Base3Lett)," (LSM Rank = ",Length(CountLSM), ").\n");
    #LSMMat := List(TransposedMat(Inverse(Mat*Z(2))),x->GF2ToZ(x));
    LSMMat := List(TransposedMat(InverseMatMod(Mat,2)));
    LSMLett := List([1..Length(CountLSM)],x->LSMMat[x]);
    Print("LSM:\n");
    #Print(Mat);
    j := 1;
    for i in [1..Length(IWP[IT])] do
        Print(IWP[IT][i][1]," ");
        if (IWP[IT][i][2] = []) = false then
            PrintMonomialString(IndToElem(LSMLett[j],Base3Lett),GensDim1to4,"+",GENNAMES[IT],"\n");
            j := j+1;
        else
            Print("\n");
        fi;
    od;
    #for j in [(Length(IWP[IT])+1)..Length(Mat)] do
    #    PrintMonomialString(IndToElem(LSMMat [j],Base3Lett),GensDim1to4,"+",GENNAMES[IT],"\n");
    #    Print(overcomplete_g[j],"\n");
    #od;
    
else
    Print("Full Rank NOT achieved: ", RankMatrix(Mat*Z(2)),"!=", Length(Base3Lett), "or", RankMatrix(Mat*Z(2)),"!=", Length(Mat),".\n");
    #Print(Mat*Z(2),"\n");
fi;
#Print("===========================================\n");



return true;
end;
#####################################################################
#####################################################################
