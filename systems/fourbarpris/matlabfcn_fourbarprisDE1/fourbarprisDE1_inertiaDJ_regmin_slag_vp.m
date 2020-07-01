% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((1+1)*1/2)x9]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbarprisDE1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaDJ_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_inertiaDJ_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaDJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:54
% EndTime: 2020-06-27 17:14:56
% DurationCPUTime: 0.19s
% Computational Cost: add. (521->32), mult. (301->72), div. (123->12), fcn. (6->2), ass. (0->40)
t28 = (qJ(1) + pkin(3));
t42 = -pkin(2) + t28;
t19 = (pkin(1) + t42);
t11 = 1 / t19;
t26 = 1 / (t28 ^ 2);
t43 = -pkin(2) - t28;
t17 = (pkin(1) - t43);
t7 = 1 / t17;
t51 = (t26 * t7);
t18 = (pkin(1) - t42);
t55 = 1 / t18;
t45 = (t55 * t51);
t20 = (pkin(1) + t43);
t54 = 1 / t20;
t40 = (t54 * t45);
t58 = (t11 * t40);
t12 = 1 / (t19 ^ 2);
t57 = (t12 * t54);
t29 = (qJ(1) ^ 2);
t53 = -2 * qJ(1);
t4 = pkin(1) ^ 2 - pkin(2) ^ 2 - t29 + (t53 - pkin(3)) * pkin(3);
t48 = qJD(1) * t4 ^ 2;
t6 = t28 * qJD(1);
t56 = 4 * t4 * t6 * t58;
t25 = 1 / t28;
t10 = 1 / (t18 ^ 2);
t14 = 1 / (t20 ^ 2);
t52 = (t7 * t55);
t49 = (t19 * t20);
t24 = t28 ^ 2;
t47 = 2 * t25 / t24 * t7;
t46 = t24 * t52;
t44 = t29 * t51;
t41 = (t18 * t49);
t39 = t55 * t44;
t38 = (t17 * t41);
t8 = 1 / (t17 ^ 2);
t37 = -t10 * t7 + t55 * t8;
t35 = t56 + (t12 * t40 + (-t14 * t45 + (t37 * t26 + t47 * t55) * t54) * t11) * t48;
t1 = [t35, 0, 0, (-4 * t6 * t25 + (-2 * t26 + t25 * (-t41 + (t49 + (t19 - t20) * t18) * t17) / t38) * t4 * qJD(1)) * (-t38) ^ (-0.1e1 / 0.2e1), 2 * t35 * qJ(1) - 2 * t48 * t58, t29 * t56 + (t39 * t57 + (-t14 * t39 - (t10 * t44 - (t29 * t47 + (t29 * t8 + t7 * t53) * t26) * t55) * t54) * t11) * t48, 4 * (t46 * t57 + (-t14 * t46 + (t37 * t24 - 2 * t28 * t52) * t54) * t11) * qJD(1), 0, 0;];
MMD_reg = t1;
