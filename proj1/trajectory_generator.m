function s_des = trajectory_generator(t, path, h)

persistent px
persistent py
persistent pz
persistent pvx
persistent pvy
persistent pvz

persistent sn
persistent tp


if nargin > 1 % pre-process can be done here (given waypoints)
    
    % generate trajectory base on constrained QP minimizing snap
    
    % write the Qk matrix, Jk = pk'*Qk*pk
    tp=2.6;
    Qk=zeros(8,8);
    for i=5:8
        for j=5:8
            Qk(i,j)=(tp^(i+j-7))*i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3)/(i+j-7);
        end
    end
    %     Qk
    
    % combine each segments' Qk matrix to Q(matrix for the whole
    % trajectory)
    ptn=size(path,1);
    sn=ptn-1;
    Q=zeros(8*sn,8*sn);
    for i=1:sn
        Q(8*(i-1)+1:8*i,8*(i-1)+1:8*i)=Qk;
    end
    
    H=zeros(24*sn,24*sn);
    for i=1:3
        H(1+8*sn*(i-1):8*sn*i,1+8*sn*(i-1):8*sn*i)=Q;
    end
    H

    % end point derivative mapping matrix, A_0_0,A_0_t ... A_2_0,A_2_t
    A_0_0=zeros(1,8);
    A_0_0(1)=1;
    
    A_0_t(1)=1;
    for i=2:8
        A_0_t(i)=A_0_t(i-1)*tp;
    end
    
    A_1_0=zeros(1,8);
    A_1_0(2)=1;
    
    A_1_t(1)=0;
    A_1_t(2)=1;
    for i=3:8
        A_1_t(i)=(i-1)*tp^(i-2);
    end
    
    A_2_0=zeros(1,8);
    A_2_0(3)=2;
    
    A_2_t=zeros(1,8);
    A_2_t(3)=2;
    for i=4:8
        A_2_t(i)=(i-1)*(i-2)*tp^(i-3);
    end
    
    % A_0_0
    % A_0_t
    % A_1_0
    % A_1_t
    % A_2_0
    % A_2_t
    
    % continuity constrains matrix
    Ak=zeros(3,16);
    Ak(1,1:16)=[A_0_t ,-A_0_0];
    Ak(2,1:16)=[A_1_t ,-A_1_0];
    Ak(3,1:16)=[A_2_t ,-A_2_0];
    
    % combine A_k_t to continuity constrains matrix of the whole trajectory
    A_contin=zeros(3*(sn-1),8*sn);
    for i=1:(sn-1)
        A_contin(1+3*(i-1):3+3*(i-1),1+8*(i-1):16+8*(i-1))=Ak;
    end
    A_contin
    
    % derivative constrains matrix
    A_deriv=zeros(ptn+4,8*sn);
    for i=1:sn
        A_deriv(i,1+8*(i-1):8+8*(i-1))=A_0_0;
    end
    A_deriv(ptn,1+8*(sn-1):8*sn)=A_0_t;
    A_deriv(ptn+1,1:8)=A_1_0;
    A_deriv(ptn+2,1:8)=A_2_0;
    A_deriv(ptn+3,1+8*(sn-1):8*sn)=A_1_t;
    A_deriv(ptn+4,1+8*(sn-1):8*sn)=A_2_t;
    
    A_deriv
    A_phi=[A_contin;A_deriv];
    
    [n,m]=size(A_phi);
    Aeq=zeros(3*n,3*m);
    Aeq(1:n,1:m)=A_phi;
    Aeq(1+n:2*n,1+m:2*m)=A_phi;
    Aeq(1+2*n:3*n,1+2*m:3*m)=A_phi;  
    Aeq
    
    % b for continuity and b for derivative
    b_contin=zeros(3*(sn-1),1);
    b_deriv_x=zeros(ptn+4,1);
    b_deriv_y=zeros(ptn+4,1);
    b_deriv_z=zeros(ptn+4,1);
    for i=1:ptn
        b_deriv_x(i)=path(i,1);
        b_deriv_y(i)=path(i,2);
        b_deriv_z(i)=path(i,3);
    end

    for i=1:4
        b_deriv_x(ptn+i)=0;
        b_deriv_y(ptn+i)=0;
        b_deriv_z(ptn+i)=0;
    end

    % complete b
    b_x=[b_contin;b_deriv_x];
    b_y=[b_contin;b_deriv_y];
    b_z=[b_contin;b_deriv_z];
    beq=[b_x;b_y;b_z];
    beq

    % create f 
    f=zeros(1,24*sn);
    f
    
    % create initial value, use straight lines
    x0=zeros(24*sn,1);
    x0_x=zeros(8*sn,1);
    x0_y=zeros(8*sn,1);
    x0_z=zeros(8*sn,1);

    for i=1:sn
        x0_x(1+8*(i-1))=path(i,1);
        x0_x(2+8*(i-1))=(path(i+1,1)-path(i,1))/tp;
        x0_y(1+8*(i-1))=path(i,2);
        x0_y(2+8*(i-1))=(path(i+1,2)-path(i,2))/tp;
        x0_z(1+8*(i-1))=path(i,3);
        x0_z(2+8*(i-1))=(path(i+1,3)-path(i,3))/tp;
    end
    x0=[x0_x;x0_y;x0_z];
    
    x=quadprog(H,f,[],[],Aeq,beq,[],[],x0);
    x
    
    %convert result to matrix, each row is one segment
    px=zeros(sn,8);
    py=zeros(sn,8);
    pz=zeros(sn,8);
    for i=1:sn
        px(i,1:8)=x(1+8*(i-1):8*i)';
        py(i,1:8)=x(8*sn+1+8*(i-1):8*sn+8*i)';
        pz(i,1:8)=x(16*sn+1+8*(i-1):16*sn+8*i)';
    end

    px
    py
    pz

    % get coefficient of velocity
    pvx=zeros(sn,8);
    pvy=zeros(sn,8);
    pvz=zeros(sn,8);
    for s=1:sn
        for i=1:7
            pvx(s,i)=i*px(s,i+1);
            pvy(s,i)=i*py(s,i+1);
            pvz(s,i)=i*pz(s,i+1);
        end
    end

else % output desired trajectory here (given time)

    % decide which segment lay in 
    segment=1;
    while 1
        if t<=segment*tp
            break
        else
            segment=segment+1;
        end
    end
    segment

    if segment>sn
        return
    end

    % get position
    time=ones(8,1);
    ts=t-(segment-1)*tp;
    for i=1:8 
        time(i)=ts^(i-1);
    end

    pos_x=px(segment,:)*time;
    pos_y=py(segment,:)*time;
    pos_z=pz(segment,:)*time;

    v_x=pvx(segment,:)*time;
    v_y=pvy(segment,:)*time;
    v_z=pvz(segment,:)*time;

    % set to s_des
    s_des=zeros(13,1);

    s_des(1)=pos_x;
    s_des(2)=pos_y;
    s_des(3)=pos_z;
    s_des(4)=v_x;
    s_des(5)=v_y;
    s_des(6)=v_z;

    ypr = [0.0, 0.0, 0.0];
    Rot = ypr_to_R(ypr);
    q_des = R_to_quaternion(Rot);
    s_des(7:10) = q_des;
    
end

end


