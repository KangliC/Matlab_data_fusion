function [Qret,accel] = ahrs_filter(Q, w_x, w_y, w_z, a_x, a_y, a_z, m_x, m_y, m_z, deltat, zeta, beta)
    persistent b_x ;
    persistent b_z;
    persistent w_bx;
    persistent w_by;
    persistent w_bz;
    mag_flag = 0;
    
    if isempty(b_x) || isnan(b_x) || isnan(w_bx) || isempty(w_bx)
        b_x = 1;
        b_z = 0;
        w_bx = 0;
        w_by = 0;
        w_bz = 0;
    end
    
    q1 = Q.q1;
    q2 = Q.q2;
    q3 = Q.q3;
    q4 = Q.q4;
    
    % compute flux in the earth frame
	q_1_q_2 = Q.q1 * Q.q2; % recompute axulirary variables
	q_1_q_3 = Q.q1 * Q.q3;
	q_1_q_4 = Q.q1 * Q.q4;
	q_3_q_4 = Q.q3 * Q.q4;
	q_2_q_3 = Q.q2 * Q.q3;
	q_2_q_4 = Q.q2 * Q.q4;    
    
    normv = norm([a_x a_y a_z]);
    a_x = a_x / normv;
	a_y = a_y / normv;
	a_z = a_z / normv;
    
    normv = norm([m_x m_y m_z]);
    if normv ~= 0
        mag_flag = 1;
        m_x = m_x / normv;
        m_y = m_y / normv;
        m_z = m_z / normv;

        h_x = 2*m_x * (0.5 - Q.q3 * Q.q3 - Q.q4 * Q.q4) + 2*m_y * (q_2_q_3 - q_1_q_4) + 2*m_z * (q_2_q_4 + q_1_q_3);
        h_y = 2*m_x * (q_2_q_3 + q_1_q_4) + 2*m_y * (0.5 - Q.q2 * Q.q2 - Q.q4 * Q.q4) + 2*m_z * (q_3_q_4 - q_1_q_2);
        h_z = 2*m_x * (q_2_q_4 - q_1_q_3) + 2*m_y * (q_3_q_4 + q_1_q_2) + 2*m_z * (0.5 - Q.q2 * Q.q2 - Q.q3 * Q.q3);
        % normali_ the flux vector to have only components in the x and z
        b_x = sqrt((h_x * h_x) + (h_y * h_y));
        b_z = h_z;
    else
        Qret = Q;
        return
    end
    
    f1 = 2*q2*q4 - 2*q1*q3 - a_x;
	f2 = 2*q1*q2 + 2*q3*q4 - a_y;
	f3 = 1 - 2*q2*q2 - 2*q3*q3 - a_z;
	f4 = 2*b_x*(0.5 - q3*q3 - q4*q4) + 2*b_z*(q2*q4 - q1*q3) - m_x;
	f5 = 2*b_x*(q2*q3 - q1*q4) + 2*b_z*(q1*q2 + q3*q4) - m_y;
	f6 = 2*b_x*(q1*q3 + q2*q4) + 2*b_z*(0.5 - q2*q2 - q3*q3) - m_z;
    
    Jg = [-2*q3 2*q4 -2*q1 2*q2;
           2*q2 2*q1  2*q4 2*q3;
           0   -4*q2 -4*q3   0];
    Jb = [-2*b_z*q3             2*b_z*q4           (-4*b_x*q3-2*b_z*q1) (-4*b_x*q4+2*b_z*q2);
          (-2*b_x*q4+2*b_z*q2) (2*b_x*q3+2*b_z*q1) (2*b_x*q2+2*b_z*q4)  (-2*b_x*q1+2*b_z*q3);
            2*b_x*q3           (2*b_x*q4-4*b_z*q2) (2*b_x*q1-4*b_z*q3)             2*b_x*q2];
    Jgb = [Jg;
           Jb];
    fg = [f1;f2;f3];
    fgb = [f1;f2;f3;f4;f5;f6];
    if mag_flag == 1
        fgrad = Jgb' * fgb;
    else
        fgrad = Jg' * fg;
    end
	% normalise the gradient to estimate the direction of the gyroscope error
	fgrad = fgrad / norm(fgrad);

	% compute angular estimated direction of the gyroscope error
	w_err_x = 2*q1*fgrad(2) - 2*q2*fgrad(1) - 2*q3*fgrad(4) + 2*q4*fgrad(3);
	w_err_y = 2*q1*fgrad(3) + 2*q2*fgrad(4) - 2*q3*fgrad(1) - 2*q4*fgrad(2);
	w_err_z = 2*q1*fgrad(4) - 2*q2*fgrad(3) + 2*q3*fgrad(2) - 2*q4*fgrad(1);
	% compute and remove the gyroscope bias
	w_bx = w_bx + w_err_x * deltat * zeta;
	w_by = w_by + w_err_y * deltat * zeta;
	w_bz = w_bz + w_err_z * deltat * zeta;
	w_x = w_x - w_bx;
	w_y = w_y - w_by;
	w_z = w_z - w_bz;
	% compute the quaternion rate measured by gyroscopes
	qDot_omega_1 = -0.5*(q2*w_x + q3*w_y + q4*w_z);
	qDot_omega_2 = 0.5*(q1*w_x + q3*w_z - q4*w_y);
	qDot_omega_3 = 0.5*(q1*w_y - q2*w_z + q4*w_x);
	qDot_omega_4 = 0.5*(q1*w_z + q2*w_y - q3*w_x);
	% integrate the estimated quaternion rate
	Q.q1 = Q.q1 + (qDot_omega_1 - (beta * fgrad(1))) * deltat;
	Q.q2 = Q.q2 + (qDot_omega_2 - (beta * fgrad(2))) * deltat;
	Q.q3 = Q.q3 + (qDot_omega_3 - (beta * fgrad(3))) * deltat;
	Q.q4 = Q.q4 + (qDot_omega_4 - (beta * fgrad(4))) * deltat;

    % normali_ quaternion
	normv = norm([Q.q1 Q.q2 Q.q3 Q.q4]);
	Q.q1 = Q.q1 / normv;
	Q.q2 = Q.q2 / normv;
	Q.q3 = Q.q3 / normv;
	Q.q4 = Q.q4 / normv;
    
    Qret = Q;
return
