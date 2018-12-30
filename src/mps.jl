struct MPS
  A
  oc::Int64
end


function dosvdtrunc(AA,m,thres=0.0)  # AA a matrix;  keep at most m states
    (u,d,v) = svd(AA)
    prob = dot(d,d)       # total probability
    #determine trunc dimension
    errorm=length(d)
    for nn=-length(d):-1
        n=-nn
        d[n]<thres ? errorm=n : break;    
    end    
    mm = min(m,length(d),errorm) # number of states to keep
    d = d[1:mm]           # middle matrix in vector form
    trunc = prob - dot(d,d)
    U = u[:,1:mm]
    V = v[:,1:mm]'
    (U,d,V,trunc)         # AA == U * diagm(d) * V   with error trunc
end

function dosvdleftright(AA,m,toright,thres=0.0)
    (U,d,V,trunc) = dosvdtrunc(AA,m,thres)
    if toright   #true -> svd to the right
    V = diagm(d) * V
    else
    U = U * diagm(d)
    end
    (U,V,trunc)
end

function dosvd4(AA,m,toright,thres=0.0)   # AA is ia * 2 * 2 * ib;  svd down the middle;  return two parts
    ia = size(AA,1)
    ib = size(AA,4)
    AA = reshape(AA,ia*2,2*ib)
    (U,V,trunc) = dosvdleftright(AA,m,toright,thres)
    mm = size(U,2)
    U = reshape(U,ia,2,mm)
    V = reshape(V,mm,2,ib)
    (U,V,trunc)
end


function totnormsq(ψ::MPS)
    totn=length(ψ.A)
    A1=copy(ψ.A[1])
    cA1=conj(ψ.A[1])
     @tensor T[cp,c]:=(A1)[a,b,c]*cA1[a,b,cp]
    for l=2:totn-1
        A1=ψ.A[l]
        cA1=conj(ψ.A[l])
        @tensor Tx[cp,c]:=T[ap,a]*A1[a,b,c]*cA1[ap,b,cp] 
        T=(Tx)
    end
    A1=ψ.A[totn]
    cA1=conj(ψ.A[totn])
    @tensor normsqval=scalar(T[ap,a]*A1[a,b,c]*cA1[ap,b,c])
    return normsqval
end

function mpsdot(ψ1::MPS,ψ2::MPS)
    totn=length(ψ1.A)
    A1=copy(ψ1.A[1])
    cA1=conj(ψ2.A[1])
     @tensor T[cp,c]:=(A1)[a,b,c]*cA1[a,b,cp]
    for l=2:totn-1
        A1=copy(ψ1.A[l])
        cA1=conj(ψ2.A[l])
        @tensor Tx[cp,c]:=T[ap,a]*A1[a,b,c]*cA1[ap,b,cp] 
        T=(Tx)
    end
    A1=ψ1.A[totn]
    cA1=conj(ψ2.A[totn])
    @tensor normsqval=scalar(T[ap,a]*A1[a,b,c]*cA1[ap,b,c])
    return normsqval
end

function mpsnorm(ψ::MPS)
    ocsite=ψ.oc;
    val=0.0;
    M=ψ.A[ocsite];
    cM=conj(M)
    @tensor val=scalar(M[a,b,c]*cM[a,b,c])
    return sqrt(real(val))
end

function moveocone!(Ai,Ai1,toright)
    (lbd,phyd,mbd)=size(Ai);
    (mbd,phyd,rbd)=size(Ai1);

    if toright
        MatAi=reshape(Ai,lbd*phyd,mbd);
        AiQ,AiR=qr(MatAi)
        @tensor newAi1[a,b,c]:=AiR[a,x]*Ai1[x,b,c]
        A1=(reshape(AiQ,lbd,phyd,mbd));
        Ai[:,:,:]=(A1)
        Ai1[:,:,:]=(newAi1);
    else
        MatAi1=reshape(Ai1,mbd,phyd*rbd);
        AiQ,AiR=qr(ctranspose(MatAi1))
        @tensor newAi[a,b,c]:=Ai[a,b,x]*ctranspose(AiR)[x,c]
        Ai1[:,:,:]=(reshape(ctranspose(AiQ),mbd,phyd,rbd));
        Ai[:,:,:]=(newAi);
        
    end
end    

function moveto!(ψ::MPS,newoc::Int64)  #move oc of an MPS to r
    oldoc=ψ.oc;
    dist=newoc-oldoc;
    if dist==0
        return ;
    end
    dist>0 ? toright=true : toright=false
    
    count=abs(dist)
    if toright
       for i=1:count
            Ai=ψ.A[oldoc+i-1]
            Ai1=ψ.A[oldoc+i]
            moveocone!(Ai,Ai1,toright)
            ψ.A[oldoc+i-1]=Ai;
            ψ.A[oldoc+i]=Ai1;
       end    
    else
         for i=1:count
            Ai=ψ.A[oldoc-i]
            Ai1=ψ.A[oldoc-i+1]
            moveocone!(Ai,Ai1,false)
            ψ.A[oldoc-i]=Ai;
            ψ.A[oldoc-i+1]=Ai1;
        end
    end    
    ψ.oc=newoc;
end

function expσz(ψ::MPS,toright=true)
    O=[1 0;0 -1];
    Ltot=length(ψ.A);
    exp=0.0;
    if toright
        for site=1:Ltot
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
        end
    else   
        for ssite=-Ltot:-1
        site=-ssite
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
        end
    end    
    return exp
end    

function expσx(ψ::MPS,toright=true)
    O=[0 1;1 0];
    Ltot=length(ψ.A);
    exp=0.0;
    if toright
       for site=1:Ltot
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
       end
    else   
        for ssite=-Ltot:-1
        site=-ssite
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
        end
    end    
    return exp
end 

function expσy(ψ::MPS,toright=true)
    O=[0 -1.0im;1.0im 0];
    Ltot=length(ψ.A);
    exp=0.0;
    if toright
       for site=1:Ltot
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
        end
    else   
        for ssite=-Ltot:-1
        site=-ssite
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval/Ltot
        end
    end    
    return exp
end

function expO1(ψ::MPS,O,toright=true)
    Ltot=length(ψ.A);
    exp=0.0+0.0im;
    if toright
       for site=1:Ltot
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
         exp+=oval
        end
    else   
        for ssite=-Ltot:-1
        site=-ssite
        moveto!(ψ,site);
        Asite=ψ.A[site]
        @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
        exp+=oval
        end
    end 
    return exp
end 

function expO1_j(ψ::MPS,O,site)
    Ltot=length(ψ.A);
    exp=0.0+0.0im;
    moveto!(ψ,site);
    Asite=ψ.A[site]
    @tensor oval=scalar(Asite[a,b,c]*conj(Asite[a,bp,c])*O[bp,b])
    exp+=oval
   return exp
end 

function expO2_j(ψ::MPS,O,j)
    L=length(ψ.A);
    if j!=L
    exp=0;
    moveto!(ψ,j);
    Aj=ψ.A[j]
    Aj1=ψ.A[j+1]
        @tensor oval=scalar(Aj[a,b,c]*conj(Aj[a,bp,cp])*Aj1[c,d,e]*conj(Aj1[cp,dp,e])*O[bp,dp,b,d])
    exp+=oval
        
    end
   return exp
end 


function expO2(ψ::MPS,O,toright=true)
    L=length(ψ.A);
    exp=0;
    if toright
       for j=1:L-1
        moveto!(ψ,j);
    Aj=ψ.A[j]
    Aj1=ψ.A[j+1]
        @tensor oval=scalar(Aj[a,b,c]*conj(Aj[a,bp,cp])*Aj1[c,d,e]*conj(Aj1[cp,dp,e])*O[bp,dp,b,d])
        exp+=oval
        end
    else   
        for ssite=-L+1:-1
        site=-ssite
        moveto!(ψ,j);
        Aj=ψ.A[j]
        Aj1=ψ.A[j+1]
        @tensor oval=scalar(Aj[a,b,c]*conj(Aj[a,bp,cp])*Aj1[c,d,e]*conj(Aj1[cp,dp,e])*O[bp,dp,b,d])
        exp+=oval
        end
    end 
    return exp
end 

