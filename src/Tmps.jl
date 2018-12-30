struct TMPS
  A
  oc::Int64
end


function totnormsq(ψ::TMPS)
    totn=length(ψ.A)
    A1=ψ.A[1]
    #cA1=conj(ψ.A[1])
     @tensor T[:]:=(A1)[1,2,3,-2]*conj(A1[1,2,3,-1])
    for l=2:totn-1
        A1=ψ.A[l]
        #cA1=conj(ψ.A[l])
        @tensor Tx[:]:=T[2,1]*A1[1,3,4,-2]*conj(A1[2,3,4,-1]) 
        T=(Tx)
    end
    A1=ψ.A[totn]
    #cA1=conj(ψ.A[totn])
    @tensor normsqval=scalar(T[2,1]*A1[1,3,4,5]*conj(A1[2,3,4,5]))
    return normsqval
end

function Tmpsdot(ψ2::TMPS,ψ1::TMPS)
    totn=length(ψ1.A)
    A1=ψ1.A[1]
    A2=ψ2.A[1]
     @tensor T[:]:=(A1)[1,2,3,-2]*conj(A2[1,2,3,-1])
    for l=2:totn-1
        A1=ψ1.A[l]
        A2=(ψ2.A[l])
        @tensor Tx[:]:=T[2,1]*A1[1,3,4,-2]*conj(A2[2,3,4,-1]) 
        T=(Tx)
    end
    A1=ψ1.A[totn]
    A2=ψ2.A[totn]
    @tensor normsqval=scalar(T[2,1]*A1[1,3,4,5]*conj(A2[2,3,4,5]))
    return normsqval
end

function Tmpsnorm(ψ::TMPS)
    ocsite=ψ.oc;
    val=0.0;
    M=ψ.A[ocsite];
    @tensor val=scalar(M[1,2,3,4]*conj(M[1,2,3,4]))
    return sqrt(real(val))
end

function moveocone_TMPS!(Ai,Ai1,toright)
    (lbd,phyd,auxd,mbd)=size(Ai);
    (mbd,phyd,auxd,rbd)=size(Ai1);

    if toright
        MatAi=reshape(Ai,lbd*phyd*auxd,mbd);
        AiQ,AiR=qr(MatAi)
        @tensor newAi1[a,b,d,c]:=AiR[a,x]*Ai1[x,b,d,c]
        Ai[:,:,:,:]=(reshape(AiQ,lbd,phyd,auxd,mbd));
        Ai1[:,:,:,:]=(newAi1);
    else
        MatAi1=reshape(Ai1,mbd,phyd*auxd*rbd);
        AiQ,AiR=qr(transpose(MatAi1))
        @tensor newAi[a,b,d,c]:=Ai[a,b,d,x]*transpose(AiR)[x,c]
        Ai1[:,:,:,:]=(reshape(transpose(AiQ),mbd,phyd,auxd,rbd));
        Ai[:,:,:,:]=(newAi);
        
    end
end    

function moveto!(ψ::TMPS,newoc::Int64)  #move oc of an MPS to r
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
            moveocone_TMPS!(Ai,Ai1,toright)
            ψ.A[oldoc+i-1]=Ai;
            ψ.A[oldoc+i]=Ai1;
       end    
    else
         for i=1:count
            Ai=ψ.A[oldoc-i]
            Ai1=ψ.A[oldoc-i+1]
            moveocone_TMPS!(Ai,Ai1,false)
            ψ.A[oldoc-i]=Ai;
            ψ.A[oldoc-i+1]=Ai1;
        end
    end    
    ψ.oc=newoc;
end


function Tmpstrace(ψ::TMPS)
    totn=length(ψ.A)
    A1=ψ.A[1]
    @tensor T[:]:=A1[-1,1,1,-2]
    for l=2:totn-1
        A1=ψ.A[l]
        @tensor Tx[:]:=T[-1,2]*A1[2,1,1,-2]
        T=(Tx)
    end
    A1=ψ.A[totn]
    @tensor traceval=scalar(T[3,2]*A1[2,1,1,3])
    return traceval
end 

function normalizeMPS!(ψ::TMPS)
    if ψ.oc==length(ψ.A)
        moveto!(ψ,length(ψ.A)-1)
    end
    oc=ψ.oc
    Ai=ψ.A[oc]
    Ai1=ψ.A[oc+1]
    (lbd,phyd,auxd,mbd)=size(Ai)
    (mbd,phyd,auxd,rbd)=size(Ai1)
    @tensor AiAi1[lb,phy1,aux1,phy2,aux2,rb]:=Ai[lb,phy1,aux1,md]*Ai1[md,phy2,aux2,rb]
    AA=reshape(AiAi1,lbd*phyd*auxd,phyd*auxd*rbd)
    #do svd find the schmidt value
    (u,d,v) = svd(AA)
    mm=mbd;
    d = d[1:mm]
    U = u[:,1:mm]
    V = v[:,1:mm]'
    #determine the norm
    totnorm2=dot(d,d)
    d=d/sqrt(totnorm2) #normalize the MPS
    U = U * diagm(d) #put the OC toleft
    #put the matrices back
    ψ.A[oc]=reshape(U,lbd,phyd,auxd,mbd);
    ψ.A[oc+1]=reshape(V,mbd,phyd,auxd,rbd);
    return totnorm2
end

