#(C) Graham Ellis, 2005-2006
#Modified by Chunxiao Liu, 2024

#Below contains three functions:

#ResolutionAlmostCrystalGroupNCH
#ResolutionExtensionNCH
#TwistedTensorProductNCH


LoadPackage("HAP");

#####################################################################
#####################################################################
TwistedTensorProductNCH:=function(R,S,EhomG,GmapE,NhomE,NEhomN,EltsE,Mult,InvE)
local
        DimensionR,BoundaryR,HomotopyR,
        DimensionS,BoundaryS,HomotopyS,
        Dimension,Boundary,Homotopy,
                FilteredLength, FilteredDimension,
        DimPQ,DimPQrec,
        Int2Pair, Pair2Int,
        Htpy, HtpyRecord, CompHtpy,
        Del, CompDel, DelRecord,
                srtfn,
        PseudoBoundary,
        Charact,
        AddWrds, AddLst,
                    HighestDegBndyMat,M2,row,hi,x,         #Added
        i,j,k,l,n,p,q,r,rr,grp,dl,
                ######Remaining variables are concerned
                ######with the contracting homotopy.
        Int2Vector,
        Vector2Int,
        HorizontalBoundaryGen,
        HorizontalBoundaryWord,
        HomotopyGradedGen,
        EmapN,
        HomotopyRec,
        Homtpy,
        HomotopyOfWord,
        FinalHomotopy,
        HorizontalPseudoBoundary,
        SMALL,SizeE,BoolE,MT,
#################################
AbsInt,                         #
SignInt;                        #
                                #
AbsInt:=AbsInt_HAP;             #
SignInt:=SignInt_HAP;           #
#################################
        
SMALL:=2^12;
if Order(S!.group) >SMALL  then SizeE:=Order(S!.group); fi;
if Order(R!.group) >SMALL
 then SizeE:=Order(R!.group); fi;
if not IsBound(SizeE) then
SizeE:=Order(R!.group)*Order(S!.group);
fi;

BoolE:=SizeE<=SMALL and Size(R!.group)=Size(R!.elts)
and Size(S!.group)=Size(S!.elts);
#################################
if BoolE then
MT:=[];
#for i in [1..SizeE] do
for i in [1..Length(EltsE)] do
MT[i]:=[];
#for j in [1..SizeE] do
for j in [1..Length(EltsE)] do
MT[i][j]:=Mult(i,j);
od;
od;
fi;
#################################

DimensionR:=R!.dimension;
DimensionS:=S!.dimension;
BoundaryR:=R!.boundary;
BoundaryS:=S!.boundary;
HomotopyR:=R!.homotopy;
HomotopyS:=S!.homotopy;
n:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")>0
then Charact:=EvaluateProperty(S,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")>0
then Charact:=Product(Intersection([
DivisorsInt(EvaluateProperty(R,"characteristic")),
DivisorsInt(EvaluateProperty(S,"characteristic"))
]));
fi;


if Charact=0 then AddWrds:=AddFreeWords; else
    AddWrds:=function(v,w);
    return AddFreeWordsModP(v,w,Charact);
    end;
fi;

#####################################################################
AddLst:=function(v,SM)
local x,ab;

for x in v do
if not IsBound(SM[x[2]]) then SM[x[2]]:=[]; fi;
ab:=AbsInt(x[1]);
if not IsBound(SM[x[2]][ab]) then
        SM[x[2]][ab]:=[ab,x[2]];
else
Unbind(  SM[x[2]][ab]);
fi;
od;

end;
#####################################################################


#####################################################################
Dimension:=function(i)
local D,j;

if i=0 then return 1; fi;

D:=0;

for j in [0..i] do
D:=D+DimensionR(j)*DimensionS(i-j);
od;

return D;


end;
#####################################################################

DimPQrec:=List([1..n+1],i->[]);
#####################################################################
DimPQ:=function(p,q)
local D,j;

if (p<0) or (q<0) then return 0; fi;
if not IsBound(DimPQrec[p+1][q+1]) then
D:=0;
for j in [0..q] do
D:=D+DimensionR(p+q-j)*DimensionS(j);
od;
DimPQrec[p+1][q+1]:=D;
fi;
return DimPQrec[p+1][q+1];
end;
#####################################################################

#####################################################################
Int2Pair:=function(i,p,q)       #Assume that x<=DimR(p)*DimS(q).
local s,r,x;
                               #The idea is that the generator f_i in F
                #corresponds to a tensor (e_r x e_s)
x:=AbsInt(i)-DimPQ(p+1,q-1);     #with e_r in R_p, e_s in S_q. If we
s:= x mod DimensionS(q);                #input i we get output [r,s].
r:=(x-s)/DimensionS(q);

if s=0 then return [SignInt(i)*r,DimensionS(q)];
else return [SignInt(i)*(r+1),s]; fi;

end;
#####################################################################

#####################################################################
Pair2Int:=function(x,p,q)
local y;                        #Pair2Int is the inverse of Int2Pair.

y:=[AbsInt(x[1]),AbsInt(x[2])];
return
SignInt(x[1])*SignInt(x[2])*((y[1]-1)*DimensionS(q)+y[2]+DimPQ(p+1,q-1));


end;
#####################################################################

HtpyRecord:=[];
for p in [0..n] do
HtpyRecord[p+1]:=[];
for q in [0..n-p] do
HtpyRecord[p+1][q+1]:=[];
for j in [1..DimPQ(p,q)] do
HtpyRecord[p+1][q+1][j]:=[];
od;
od;
od;

if not BoolE then
#####################################################################
Htpy:=function(p,q,x)
local tensor, t,g, r, s,AB;

AB:=AbsInt(x[1]);

if not IsBound(HtpyRecord[p+1][q+1][AB][x[2]]) then

tensor:=Int2Pair(AB,p,q);


g:=NEhomN(Mult(InvE(GmapE(EhomG(x[2]))),x[2] ));
t:=GmapE(EhomG(x[2]));
r:=ShallowCopy(HomotopyS(q,[tensor[2],g]));

Apply(r,y->[y[1],NhomE(y[2])]);
Apply(r,y->[Pair2Int([tensor[1],y[1]],p,q+1),Mult(t,y[2])]);

HtpyRecord[p+1][q+1][AB][x[2]]:=r;
fi;

if SignInt(x[1])>0 then
return HtpyRecord[p+1][q+1][AB][x[2]];
else
return NegateWord(HtpyRecord[p+1][q+1][AB][x[2]]);
fi;

end;
#####################################################################
else
#####################################################################
Htpy:=function(p,q,x)
local tensor, t,g, r, s,AB;

AB:=AbsInt(x[1]);

if not IsBound(HtpyRecord[p+1][q+1][AB][x[2]]) then

tensor:=Int2Pair(AB,p,q);

g:=NEhomN(MT[InvE(GmapE(EhomG(x[2])))][x[2]] );
t:=GmapE(EhomG(x[2]));
r:=ShallowCopy(HomotopyS(q,[tensor[2],g]));
Apply(r,y->[y[1],NhomE(y[2])]);
Apply(r,y->[Pair2Int([tensor[1],y[1]],p,q+1),MT[t][y[2]]]);
HtpyRecord[p+1][q+1][AB][x[2]]:=r;
fi;

if SignInt(x[1])>0 then
return HtpyRecord[p+1][q+1][AB][x[2]];
else
return NegateWord(HtpyRecord[p+1][q+1][AB][x[2]]);
fi;

end;
#####################################################################
fi;




if not Charact=2 then
#####################################################################
CompHtpy:=function(p,q,b)
local w, r;

r:=[];
for w in b do
r:=AddWrds(Htpy(p,q,w),r);
od;

return r;
end;
#####################################################################
else
#####################################################################
CompHtpy:=function(p,q,b)
local w, x, r, SM,j;

SM:=[];
for w in b do
AddLst(Htpy(p,q,w),SM);
od;

r:=[];
for x in SM do
for j in x do
if not j=0 then
Add(r,j);fi;
od;
od;

return r;
end;
#####################################################################
fi;
#####################################################################
Del:=function(k,p,q,x)
local b,i,r,v,w,tensor,Record,Ab,j,SM,y;  #Assume that 1 <= x <= DimR(p)*DimS(q)

Ab:=AbsInt(x);
if not DelRecord[k+1][p+1][q+1][Ab] = 0 then
    if SignInt(x)=1 then return DelRecord[k+1][p+1][q+1][Ab];
    else return NegateWord( DelRecord[k+1][p+1][q+1][Ab] );
    fi;
fi;

    #############################################################
    Record:=function();
    if SignInt(x)=1 then
    DelRecord[k+1][p+1][q+1][Ab]:=v;
    else
    DelRecord[k+1][p+1][q+1][Ab]:=NegateWord(v);
    fi;
    end;
    #############################################################

tensor:=Int2Pair(x,p,q);

if k=0 then
   b:=BoundaryS(q,tensor[2]);
   v:= List(b,v->[Pair2Int([tensor[1],v[1]],p,q-1),NhomE(v[2])]);
   Record();
   return v;
fi;

if k=1 then
   if q=0 then
         if p>0 then
            b:=ShallowCopy(BoundaryR(p,-tensor[1]));
            v:=List(b,v->[Pair2Int([v[1],tensor[2]],p-1,q),GmapE(v[2])]);
        Record();
        return v;
     else return [];
     fi;
   else
         if p>0 then
            v:=CompHtpy(p-1,q-1,CompDel(1,p,q-1,Del(0,p,q,-x)));
        Record();
        return v;
         else return [];
     fi;
   fi;
fi;

if k>1 then
    if p>(k-1) then

       r:=[];
       for i in [1..k] do
###
### THE NEXT LINE TAKES UP ALL THE TIME!
       r:=AddWrds(CompDel(i,p-k+i,q+k-i-1,Del(k-i,p,q,-x)),r);
###
###
       od;
       v:= CompHtpy(p-k,q+k-2,r);
       Record();
       return v;
    else return [];
    fi;
fi;

end;
#####################################################################


if not BoolE then
#####################################################################
CompDel:=function(k,p,q,b)
local r,v,w,x, map,SM,j,y;

map:=function(x);
return Del(k,p,q,x);
end;

###############
if not Charact=2 then
r:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],Mult(v[2],y[2])]);
r:=AddWrds(w,r);
od;

else
SM:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],Mult(v[2],y[2])]);
AddLst(w,SM);
od;

r:=[];
for y in SM do
for j in y do
#if not j=0 then
Add(r,j);
#fi;
od;
od;
fi;
###############

return r;
end;
#####################################################################
else
#####################################################################
CompDel:=function(k,p,q,b)
local r,v,w,x, map,SM,j,y;

map:=function(x);
return Del(k,p,q,x);
end;

###############
if not Charact=2 then
r:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],MT[v[2]][y[2]]]);
r:=AddWrds(w,r);
od;

else
SM:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],MT[v[2]][y[2]]]);
AddLst(w,SM);
od;

r:=[];
for y in SM do
for j in y do
#if not j=0 then
Add(r,j);
#fi;
od;
od;
fi;
###############

return r;
end;
#####################################################################
fi;

DelRecord:=[];
for l in [0..n] do
DelRecord[l+1]:=[];
for p in [0..n] do
DelRecord[l+1][p+1]:=[];
for q in [0..n-p] do
DelRecord[l+1][p+1][q+1]:=[];
for j in [DimPQ(p+1,q-1)+1..DimPQ(p,q)] do
DelRecord[l+1][p+1][q+1][j]:=0;
od;
od;
od;
od;



#if n <= 3 then

PseudoBoundary:=[];
HorizontalPseudoBoundary:=[];
#for k in [1..n] do
#PseudoBoundary[k]:=[];
#HorizontalPseudoBoundary[k]:=[];
#   for q in [0..k] do
#   p:=k-q;
#      for j in [DimPQ(p+1,q-1)+1..DimPQ(p,q)] do
#      r:=[];rr:=[];
#         for l in [0..p] do
#     dl:=Del(l,p,q,j);
#         r:=AddWrds(dl,r);
#     if l>0 then rr:=AddWrds(dl,rr);fi;
#         od;
#         Add(PseudoBoundary[k],[r,p]);    #I'm pretty sure it's p and not q
#     Add(HorizontalPseudoBoundary[k],rr);
#      od;
#  od;
#od;
#fi;

#####################################################################
Boundary:=function(k,j);
if k=0 then return [];
else
    if SignInt(j)=1 then return PseudoBoundary[k][j][1];
    else return NegateWord(PseudoBoundary[k][-j][1]);
    fi;
fi;
end;
#####################################################################




hi:=6;
M2:=[];


for q in [3] do

    p:=hi+1-q;
        
    for j in [DimPQ(p+1,q-1)+1..DimPQ(p,q)] do

        row:=List([1..Dimension(hi)],x->0);
        r:=[];

        for l in [0..p] do
            dl:=Del(l,p,q,j);
            r:=AddWrds(dl,r);
        od;

        for x in r do
            row[AbsoluteValue(x[1])]:= RemInt(row[AbsoluteValue(x[1])] + 1, 2);
        od;
        #Print(row,"\n");
        Add(M2,row);
    od;
od;

Print(M2,"\n");



########START WORKING ON THE CONTRACTING HOMOTOPY####################
FinalHomotopy:=fail;
##########FINISHED WORKING ON THE CONTRACTING HOMOTOPY##############


grp:=Group(EltsE);



################spectral sequence requirements##################

FilteredLength:=Length(R);

##################################################
FilteredDimension:=function(r,i);
return Length(Filtered(List(PseudoBoundary[i],x->x[2]),y->y<=r));
end;
##################################################


return Objectify(HapResolution,
        rec(
        dimension:=Dimension,
            filteredDimension:=FilteredDimension,
        boundary:=Boundary,
        homotopy:=FinalHomotopy,
        elts:=EltsE,
        group:=grp,
        #vectorToInt:=Vector2Int,
        #intToVector:=Int2Vector,
            #pseudoBoundary:=PseudoBoundary, #modified by Chunxiao Liu -- just to save some memory
            HighestDegBndyMat:=M2, #this is the highest degree matrix!
        properties:=
        [["type","resolution"],
         ["length",n],
             ["filtration_length",FilteredLength],
             ["initial_inclusion",false],
         ["characteristic",Charact],
         ["isTwistedTensorProduct",true]
          ]));
end;
#####################################################################
#####################################################################









#####################################################################
#####################################################################
ResolutionExtensionNCH:=function(arg)
local
    EEhomGG, RN, RG, TestFinite,PreImRep,
    N,E,G,
    NhomE,
    EhomG,
    GmapE,
    NEhomN,
    NEhomNrecord,
    EltsE,
    MultE,
    InvE,
    PreimagesRecordG,PreimagesRecordE,
    NisFinite,GisFinite,EisFinite,
    Lngth,T,
    AppendToElts,
    gn,i,j,x,y;

T:=0;

EEhomGG:=arg[1];
RN:=arg[2];
RG:=arg[3];
TestFinite:=false;
if Length(arg)>3 then
if arg[4]="TestFiniteness" then TestFinite:=true; fi;
fi;
if Length(arg)>4 then
PreImRep:=arg[5];
else
    #############################################################
    PreImRep:=function(x);
    return PreImagesRepresentative(EEhomGG,x);
    end;
    #############################################################
fi;

N:=RN!.group;
E:=Source(EEhomGG);
G:=Image(EEhomGG);

NisFinite:=false;
GisFinite:=false;
EisFinite:=false;
if TestFinite then
if IsFinite(N) then
    if Order(N)<=Length(RN!.elts) then NisFinite:=true; fi;
fi;

if IsFinite(G) then
        if Order(G)<=Length(RG!.elts) then GisFinite:=true; fi;
fi;

EisFinite:=IsFinite(E);
fi;

if EisFinite then EltsE:=Elements(E);
else
EltsE:=[Identity(E)];
for gn in GeneratorsOfGroup(E) do
Append(EltsE,[gn,gn^-1]);
od;
fi;

EltsE:=SSortedList(EltsE);
gn:=Position(EltsE,Identity(E));
EltsE[gn]:=EltsE[1];
EltsE[1]:=Identity(E);


    ########################################################
    AppendToElts:=function(x);
    Append(EltsE,[x]);
    end;
    ########################################################

if GisFinite then
    #########################################
    EhomG:=function(x);
    return Position(RG!.elts,ImageElm(EEhomGG,EltsE[x]));
    end;
    #########################################
else
        #########################################
        EhomG:=function(x)
    local g,Eltg;
    Eltg:=ImageElm(EEhomGG,EltsE[x]);
        g:=Position(RG!.elts,Eltg);
    if g=fail then
    RG!.appendToElts(Eltg);
    Append(RG!.elts,[Eltg]);
    g:=Length(RG!.elts); fi;
    #if Position(RG!.elts,Eltg^-1)=fail then
    #Append(RG!.elts,[Eltg^-1]);fi;
    return g;
        end;
        #########################################
fi;

if EisFinite then
    #########################################
    NhomE:=function(x);
    return Position(EltsE,RN!.elts[x]);
    end;
    #########################################
else
        #########################################
        NhomE:=function(x)
    local e,Elte;
    Elte:=RN!.elts[x];
        e:=Position(EltsE,Elte);
    if e=fail then AppendToElts(Elte);
    e:=Length(EltsE); fi;
    #if Position(EltsE,Elte^-1)=fail then
    #Append(EltsE,[Elte^-1]); fi;
    return e;
        end;
        #########################################
fi;

if GisFinite and EisFinite then
PreimagesRecordE:=List([1..Order(G)],x->
    Position(EltsE,PreImRep(RG!.elts[x])));
    
    #########################################
    GmapE:=function(x);
    return PreimagesRecordE[x];
    end;
    #########################################
else
PreimagesRecordG:=[];
PreimagesRecordE:=[];
    #########################################
    GmapE:=function(x)
    local e,Elte,Eltg,pos;
    Eltg:=RG!.elts[x];
    pos:=Position(PreimagesRecordG,Eltg);
    if not pos=fail then
    return PreimagesRecordE[pos]; fi;
    
    Elte:=PreImRep(Eltg);
    e:=Position(EltsE,Elte);
    if e=fail then AppendToElts(Elte);
        e:=Length(EltsE); fi;
    if Position(EltsE,Elte^-1)=fail then
        AppendToElts(Elte^-1); fi;
    Append(PreimagesRecordG,[Eltg]);
    Append(PreimagesRecordE,[e]);
        return e;
    end;
    #########################################
fi;

if NisFinite then
    #########################################
    NEhomN:=function(x);
    return Position(RN!.elts,EltsE[x]);
    end;
    #########################################
else
    #########################################
    NEhomN:=function(x)
    local p,Eltp;
    Eltp:= EltsE[x];
    p:=Position(RN!.elts,Eltp);
    if p=fail then RN!.appendToElts(Eltp);
    Append(RN!.elts,[Eltp]);
    p:=Length(RN!.elts); fi;
    return p;
    end;
    #########################################
fi;

if EisFinite then
    #########################################
    MultE:=function(x,y);
    return Position(EltsE,EltsE[x]*EltsE[y]);
    end;
    #########################################
else
    #########################################
    MultE:=function(x,y)
    local p,Eltp;
    Eltp:=EltsE[x]*EltsE[y];
    p:= Position(EltsE,Eltp);
    if p=fail then AppendToElts(Eltp);
    p:=Length(EltsE); fi;
    #if Position(EltsE,Eltp^-1)=fail then
        #Append(EltsE,[Eltp^-1]); fi;
    return p;
    end;
    #########################################
fi;

if EisFinite then
    #########################################
    InvE:=function(x);
    return Position(EltsE,EltsE[x]^-1);
    end;
    #########################################
else
    #########################################
    InvE:=function(x)
    local p;
    p:=(Position(EltsE,EltsE[x]^-1));
    if p=fail then AppendToElts(EltsE[x]^-1);
    p:=Length(EltsE);fi;
    return p;
    end;
    #########################################
fi;

if (not EisFinite ) and (not Length(RN!.elts)=infinity) and HAPconstant<50 then

for x in RN!.elts do
for y in RG!.elts do
AppendToElts(x*PreImRep(y));
od;
od;

fi;

#Print("\n",[NisFinite, GisFinite, EisFinite],"\n");


T:=TwistedTensorProductNCH(RG,RN,EhomG,GmapE,NhomE,NEhomN,EltsE,MultE,InvE);

        ########################################################
        AppendToElts:=function(x);
        Append(T!.elts,[x]);
        end;
        ########################################################


T!.appendToElts:=AppendToElts;

return T;
end;
#####################################################################

#####################################################################






#####################################################################
#####################################################################
ResolutionAlmostCrystalGroupNCH:=function(G,K)
local
    GhomP,P,T,RP,RT,RG,PC,PhomPC;

if not (IsAlmostCrystallographic(G) and IsPcpGroup(G)) then
Print("This function can only be applied to Almost Crystallographic pcp groups. \n"); return fail;
fi;

GhomP:=NaturalHomomorphismOnHolonomyGroup(G);
P:=Image(GhomP);
T:=Kernel(GhomP);

PhomPC:=IsomorphismPcGroup(P);
PC:=Image(PhomPC); #For some reason, good resolutions of finite pc groups
           #seem to be easier to construct than for finite pcp groups.

RP:=ResolutionFiniteGroup(PC,K);

RP!.group:=P;
RP!.elts:=List(RP!.elts,x->PreImageElm(PhomPC,x));
RT:=ResolutionNilpotentGroup(T,K);

return ResolutionExtensionNCH(GhomP,RT,RP,"Don't Test Finiteness");
end;
#####################################################################
#####################################################################

T1228:=[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
T2228:=[[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]];
T3228:=[[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]];
C2228:=[[0,1,0,0],[1,0,0,0],[-1,-1,-1,1/2],[0,0,0,1]];
C2p228:=[[0,0,1,0],[-1,-1,-1,1/2],[1,0,0,0],[0,0,0,1]];
C3228:=[[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]];
M228:=[[0,1,0,1/2],[1,0,0,1/2],[0,0,1,1/2],[0,0,0,1]];
P228:=[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]];
G228:=Group(T1228,T2228,T3228,C2228,C2p228,C3228,M228,P228);
Gp228:=IsomorphismPcpGroup(AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G228))));
R228:=ResolutionAlmostCrystalGroupNCH(Image(Gp228),7);

