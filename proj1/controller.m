function [F, M] = controller(t, s, s_des)

global params

m = params.mass;
g = params.grav;
I = params.I;

F = 1.0; M = [0.0, 0.0, 0.0]'; % You should calculate the output F and M

persistent init;
init = false;
persistent time;
time = 0.0;
persistent last_p_d_des;
persistent last_phi_c;
persistent last_theta_c;
persistent last_psi_des;

% 1,2,3 : x y z; 4 5 6: xdot ydot zdot; 7 8 9 10: qw qx qy qz; 11 12 13: p q r(wx wy wz)

% time step
dt = t - time;
time = t;

% calculate p_dd_c
kd1 = [16.0,16.0,16.0]';
kp1 = [10.0,10.0,10.0]';

p_d_des = s_des(4:6);
p_d = s(4:6);
p_des = s_des(1:3);
p = s(1:3);

p_dd_des = [0.0;0.0;0.0];
if ~init
    p_dd_des = [0.0;0.0;0.0];
else
    p_dd_des = (p_d_des - last_p_d_des)/dt;
end
last_p_d_des = p_d_des;

p_dd_c = p_dd_des + kd1.*(p_d_des-p_d)+kp1.*(p_des-p);
disp('----------------------------------------');
disp('p:');
p'
% p_des
% p_d
% p_d_des
% p_dd_des
p_dd_c'

% calculate F
F = m*(g+p_dd_c(3));

% calculate yaw
% [phi,theta,psi] = RotToRPY_ZXY(quaternion_to_R(s(7:10))');
% [phi_des,theta_des,psi_des] = RotToRPY_ZXY(quaternion_to_R(s_des(7:10))');
[phi,theta,psi] = RotToRPY_ZXY(QuatToRot(s(7:10)));
[phi_des,theta_des,psi_des] = RotToRPY_ZXY(QuatToRot(s_des(7:10)));
% disp('angle:');
% phi
% theta
% psi
% phi_des
% theta_des
% psi_des

% calculate phi_c and theta_c
phi_c = 1/g*(p_dd_c(1)*sin(psi)-p_dd_c(2)*cos(psi));
theta_c = 1/g*(p_dd_c(1)*cos(psi)+p_dd_c(2)*sin(psi));
disp('angle c:')
phi_c'
theta_c'

phi_d_c=0.0;
theta_d_c=0.0;
psi_d_c=0.0;
if ~init 
    phi_d_c=0.0;
    theta_d_c=0.0;
    psi_d_c=0.0;
    init=true;
else
    phi_d_c=(phi_c-last_phi_c)/dt;
    theta_d_c=(theta_c-last_theta_c)/dt;
    psi_d_c=(psi_des-last_psi_des)/dt;
end
last_phi_c=phi_c;
last_theta_c=theta_c;
last_psi_des=psi_des;
% disp('angle d c:');
% phi_d_c
% theta_d_c
% psi_d_c

omega=s(11:13);
phi_d=omega(1);
theta_d=omega(2);
psi_d=omega(3);
% disp('omega:');
% omega

% calculate phi_dd_c, theta_dd_c, psi_dd_c
phi_dd_c=600*(phi_c-phi)+40.0*(phi_d_c-phi_d);
theta_dd_c=600*(theta_c-theta)+40.0*(theta_d_c-theta_d);
psi_dd_c=600*(psi_des-psi)+40.0*(psi_d_c-psi_d);
% disp('angle dd c:');
% phi_dd_c
% theta_dd_c
% psi_dd_c

% calculate M
M=I*[phi_dd_c;theta_dd_c;psi_dd_c]+cross(omega,I*omega);
disp('M:');
M'

end
















