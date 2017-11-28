#functions for doing sweeping---------------------------------------------------------
function odd_sweep!(ψ::TMPS,gateL,gateM,gateR,maxm,thres=0.0)  #implicitly always sweep to the right
    #cgateL=conj(gateL);
    #cgateM=conj(gateM);
    #cgateR=conj(gateR);
    #apply gate (HL) to MPS1 MPS2
    moveto!(ψ,1) #move the oc to 1st site to start
    Ai=ψ.A[1]
    Ai1=ψ.A[2]
    (lbd,phyd,auxd,mbd)=size(Ai);
    (mbd,phyd,auxd,rbd)=size(Ai1);
    @tensor AA[:]:=gateL[-2,-4,2,3]*Ai[-1,2,4,1]*Ai1[1,3,5,-6]*conj(gateL[4,5,-3,-5]);
    (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),maxm[1],true,thres)
    (trash,mbd)=size(U)
    ψ.A[1]=reshape(U,lbd,phyd,auxd,mbd)
    ψ.A[2]=reshape(V,mbd,phyd,auxd,rbd)
    ψ.oc=2;
    #apply gate (HM) in the middle
    for i=3:2:L-2
        moveto!(ψ,i)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
    (lbd,phyd,auxd,mbd)=size(Ai);
    (mbd,phyd,auxd,rbd)=size(Ai1);
      @tensor AA[:]:=gateM[-2,-4,2,3]*Ai[-1,2,4,1]*Ai1[1,3,5,-6]*conj(gateM[4,5,-3,-5]);
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),maxm[i],true,thres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i+1;
    end
    #apply gate(HR) in the end
    i=L-1;   
        moveto!(ψ,i)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
       (lbd,phyd,auxd,mbd)=size(Ai);
       (mbd,phyd,auxd,rbd)=size(Ai1);
        @tensor AA[:]:=gateR[-2,-4,2,3]*Ai[-1,2,4,1]*Ai1[1,3,5,-6]*conj(gateR[4,5,-3,-5]);
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),maxm[i],true,thres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i+1;
end

function even_sweep!(ψ::TMPS,gateM2,maxm,thres=0.0)  #implicitly always sweep to the left
    for ii=-(L-1)+1:2:-2
        i=-ii;
        moveto!(ψ,i+1)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
       (lbd,phyd,auxd,mbd)=size(Ai);
       (mbd,phyd,auxd,rbd)=size(Ai1);
        @tensor AA[:]:=gateM2[-2,-4,2,3]*Ai[-1,2,4,1]*Ai1[1,3,5,-6]*conj(gateM2[4,5,-3,-5]);
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),maxm[i],false,thres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i;
    end
end

function imag_odd_sweep!(ψ::TMPS,gateL,gateM,gateR,βmaxm,βthres=0.0)  #implicitly always sweep to the right

    #apply gate (HL) to MPS1 MPS2
    moveto!(ψ,1) #move the oc to 1st site to start
    Ai=ψ.A[1]
    Ai1=ψ.A[2]
    (lbd,phyd,auxd,mbd)=size(Ai);
    (mbd,phyd,auxd,rbd)=size(Ai1);
    @tensor AA[:]:=gateL[-2,-4,2,3]*Ai[-1,2,-3,1]*Ai1[1,3,-5,-6];
    (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),βmaxm[1],true,βthres)
    (trash,mbd)=size(U)
    ψ.A[1]=reshape(U,lbd,phyd,auxd,mbd)
    ψ.A[2]=reshape(V,mbd,phyd,auxd,rbd)
    ψ.oc=2;
    #apply gate (HM) in the middle
    for i=3:2:L-2
        moveto!(ψ,i)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
    (lbd,phyd,auxd,mbd)=size(Ai);
    (mbd,phyd,auxd,rbd)=size(Ai1);
        @tensor AA[:]:=gateM[-2,-4,2,3]*Ai[-1,2,-3,1]*Ai1[1,3,-5,-6];
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),βmaxm[i],true,βthres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i+1;
    end
    #apply gate(HR) in the end
    i=L-1;   
        moveto!(ψ,i)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
       (lbd,phyd,auxd,mbd)=size(Ai);
       (mbd,phyd,auxd,rbd)=size(Ai1);
        @tensor AA[:]:=gateR[-2,-4,2,3]*Ai[-1,2,-3,1]*Ai1[1,3,-5,-6];
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),βmaxm[i],true,βthres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i+1;
end

function imag_even_sweep!(ψ::TMPS,gateM2,βmaxm,βthres=0.0)  #implicitly always sweep to the left
   
    for ii=-(L-1)+1:2:-2
        i=-ii;
        moveto!(ψ,i+1)
        Ai=ψ.A[i]
        Ai1=ψ.A[i+1]
       (lbd,phyd,auxd,mbd)=size(Ai);
       (mbd,phyd,auxd,rbd)=size(Ai1);
        @tensor AA[a,p1p,q1,p2p,q2,b]:=gateM2[-2,-4,2,3]*Ai[-1,2,-3,1]*Ai1[1,3,-5,-6];
        (U,V,trunc)=dosvdleftright(reshape(AA,lbd*phyd*auxd,phyd*auxd*rbd),βmaxm[i],false,βthres)
        (trash,mbd)=size(U)
        ψ.A[i]=reshape(U,lbd,phyd,auxd,mbd)
        ψ.A[i+1]=reshape(V,mbd,phyd,auxd,rbd)
        ψ.oc=i;
    end
end