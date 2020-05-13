% Jacobian of implicit kinematic constraints of
% mg10hlIC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:09
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = mg10hlIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:09:18
% EndTime: 2020-04-11 13:09:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (230->30), mult. (262->46), div. (14->2), fcn. (196->12), ass. (0->75)
unknown=NaN(7,6);
t1 = cos(qJ(10));
t2 = t1 * pkin(4);
t3 = qJ(2) + pkin(14) + qJ(3);
t4 = cos(t3);
t5 = t4 * pkin(3);
t7 = sin(t3);
t8 = t7 * pkin(3);
t9 = sin(qJ(10));
t10 = t9 * pkin(4);
t13 = 0.1e1 / (t8 * t10 + t2 * t5);
t14 = qJ(2) + pkin(14);
t15 = cos(t14);
t17 = pkin(1) * t15 + t5;
t18 = t13 * t17;
t20 = sin(t14);
t22 = pkin(1) * t20 + t8;
t26 = cos(qJ(5));
t27 = t26 * pkin(7);
t28 = qJ(9) + qJ(11);
t29 = cos(t28);
t30 = t29 * qJ(13);
t32 = cos(qJ(9));
t33 = t32 * pkin(5);
t34 = t30 - t33;
t35 = sin(qJ(5));
t38 = t35 * pkin(7);
t40 = sin(t28);
t41 = t40 * qJ(13);
t42 = sin(qJ(9));
t43 = t42 * pkin(5);
t44 = -t41 + t43;
t51 = 0.1e1 / (-t34 * t35 * pkin(7) - t44 * t29 * qJ(13) - t27 * t41 - t27 * t44 + t30 * t38 - t41 * t34);
t57 = -(t27 + t30) * t51 * t29 + (-t38 - t41) * t51 * t40;
unknown(1,1) = 0;
unknown(1,2) = (-t10 * t13 * t22 - t2 * t18);
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(2,1) = 0;
unknown(2,2) = 0;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = 0;
unknown(2,6) = t57;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = 0;
unknown(3,6) = (t33 * t51 * t29 + t43 * t51 * t40);
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 0;
unknown(4,5) = 0;
unknown(4,6) = -t57;
unknown(5,1) = 0;
unknown(5,2) = (-t13 * t4 * pkin(3) * t22 + t8 * t18);
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(5,5) = 0;
unknown(5,6) = 0;
unknown(6,1) = 0;
unknown(6,2) = 0;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(6,5) = 0;
unknown(6,6) = (-(t27 + t30 - t33) * t51 * t29 + (-t38 - t41 + t43) * t51 * t40);
unknown(7,1) = 0;
unknown(7,2) = (-(t8 + t2) * t13 * t17 + t13 * (t5 - t10) * t22 + 0.1e1);
unknown(7,3) = 0;
unknown(7,4) = 0;
unknown(7,5) = 0;
unknown(7,6) = 0;
B21  = unknown;
