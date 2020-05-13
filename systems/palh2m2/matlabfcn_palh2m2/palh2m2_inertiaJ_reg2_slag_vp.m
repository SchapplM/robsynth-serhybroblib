% Calculate inertial parameters regressor of joint inertia matrix for
% palh2m2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t31 = cos(qJ(2));
t47 = 0.2e1 * t31;
t28 = sin(qJ(2));
t46 = pkin(4) * t28;
t27 = sin(qJ(3));
t45 = pkin(5) * t27;
t30 = cos(qJ(3));
t44 = pkin(5) * t30;
t43 = t31 * pkin(4);
t22 = t28 ^ 2;
t25 = t31 ^ 2;
t16 = t22 + t25;
t21 = t27 ^ 2;
t24 = t30 ^ 2;
t10 = (t21 + t24) * t16;
t26 = sin(qJ(4));
t42 = t10 * t26;
t29 = cos(qJ(4));
t41 = t10 * t29;
t32 = pkin(5) ^ 2;
t40 = t21 * t32;
t39 = t27 * t16;
t38 = t30 * t16;
t37 = t26 ^ 2 + t29 ^ 2;
t12 = (t27 * t28 + t30 * t31) * pkin(4);
t13 = -t27 * t43 + t30 * t46;
t7 = -t27 * t12 - t30 * t13;
t36 = t7 * t45;
t35 = t10 * t45;
t34 = t29 * t45;
t18 = -pkin(1) - t43;
t8 = (-pkin(2) - t44) * t16 + t18;
t19 = t24 * t32;
t15 = t16 ^ 2;
t11 = -t16 * pkin(2) + t18;
t9 = t10 ^ 2;
t6 = t30 * t12 - t27 * t13;
t5 = t7 ^ 2;
t4 = t6 ^ 2;
t3 = t6 * t44;
t2 = -t10 * pkin(3) + t8;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t28 * t47, 0, t25, 0, 0, pkin(1) * t47, -0.2e1 * pkin(1) * t28, 0, pkin(1) ^ 2, 0, 0, 0, t15, 0, 0, -0.2e1 * t18 * t16, 0, 0, t18 ^ 2, t21 * t15, 0.2e1 * t27 * t15 * t30, 0, t24 * t15, 0, 0, -0.2e1 * t11 * t38, 0.2e1 * t11 * t39, 0, t11 ^ 2, 0, 0, 0, t9, 0, 0, -0.2e1 * t8 * t10, 0, 0, t8 ^ 2, 0, 0, 0, 0, 0, t9, -0.2e1 * t2 * t41, 0.2e1 * t2 * t42, 0, t37 * t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t10, 0, 0, 0, 0, 0, 0, 0, -t7 * t42, -t7 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * pkin(4) ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t13 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 + t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t5 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, t26 * t35, t10 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 - t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t36 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t40 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t29 * t2, t26 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t7, -t29 * t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t45, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t1;
