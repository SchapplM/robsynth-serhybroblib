% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% MMD_reg [((1+1)*1/2)x(1*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbarprisDE1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaDJ_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_inertiaDJ_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaDJ_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:10:11
% EndTime: 2020-05-07 09:10:14
% DurationCPUTime: 0.23s
% Computational Cost: add. (629->33), mult. (365->74), div. (153->12), fcn. (6->2), ass. (0->38)
t29 = (qJ(1) + pkin(3));
t41 = -pkin(2) + t29;
t20 = (pkin(1) + t41);
t12 = 1 / t20;
t27 = 1 / (t29 ^ 2);
t42 = -pkin(2) - t29;
t18 = (pkin(1) - t42);
t8 = 1 / t18;
t50 = (t27 * t8);
t21 = (pkin(1) + t42);
t53 = 1 / t21;
t44 = (t53 * t50);
t39 = (t12 * t44);
t19 = (pkin(1) - t41);
t54 = 1 / t19;
t56 = (t39 * t54);
t30 = (qJ(1) ^ 2);
t52 = -2 * qJ(1);
t5 = pkin(1) ^ 2 - pkin(2) ^ 2 - t30 + (t52 - pkin(3)) * pkin(3);
t47 = qJD(1) * t5 ^ 2;
t7 = t29 * qJD(1);
t55 = 4 * t5 * t7 * t56;
t26 = 1 / t29;
t11 = 1 / (t19 ^ 2);
t15 = 1 / (t21 ^ 2);
t25 = t29 ^ 2;
t51 = t25 * t8;
t48 = (t20 * t21);
t46 = 2 * t26 / t25 * t8;
t45 = t53 * t51;
t43 = t30 * t50;
t40 = (t19 * t48);
t38 = t54 * t43;
t37 = (t18 * t40);
t13 = 1 / (t20 ^ 2);
t9 = 1 / (t18 ^ 2);
t1 = t55 + (-t11 * t39 + (t13 * t44 + (-t15 * t50 + (t27 * t9 + t46) * t53) * t12) * t54) * t47;
t2 = [0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, (-4 * t7 * t26 + (-2 * t27 + t26 * (-t40 + (t48 + (t20 - t21) * t19) * t18) / t37) * t5 * qJD(1)) * (-t37) ^ (-0.1e1 / 0.2e1), 0, 2 * qJ(1) * t1 - 2 * t47 * t56, t30 * t55 + (t53 * t13 * t38 + (-t15 * t38 - (t11 * t43 - (t30 * t46 + (t30 * t9 + t52 * t8) * t27) * t54) * t53) * t12) * t47, 0, 0, 0, 0, 0, 4 * (-t11 * t12 * t45 + (t13 * t45 + (-t15 * t51 + (t25 * t9 - 2 * t29 * t8) * t53) * t12) * t54) * qJD(1), 0, 0, 0, 0;];
MMD_reg = t2;
