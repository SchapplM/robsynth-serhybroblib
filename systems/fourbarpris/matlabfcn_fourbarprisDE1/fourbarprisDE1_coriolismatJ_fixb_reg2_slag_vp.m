% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(1*1)x(1*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbarprisDE1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:10:16
% EndTime: 2020-05-07 09:10:17
% DurationCPUTime: 0.27s
% Computational Cost: add. (272->111), mult. (621->152), div. (29->7), fcn. (2->2), ass. (0->50)
t33 = (pkin(3) ^ 2);
t82 = 2 * pkin(3) * qJ(1) + t33;
t41 = pkin(3) * t33;
t44 = t41 ^ 2;
t28 = qJ(1) ^ 2;
t46 = qJ(1) * t28;
t49 = t46 ^ 2;
t76 = (pkin(1) - pkin(2));
t53 = (t76 ^ 2);
t75 = (pkin(1) + pkin(2));
t55 = (t75 ^ 2);
t81 = t76 * t53 * t75 * t55 - t44 - t49;
t20 = qJ(1) + pkin(3);
t80 = 1 / t20 ^ 2;
t37 = (pkin(1) ^ 2);
t38 = (t37 ^ 2);
t79 = 3 * t38;
t42 = t33 ^ 2;
t78 = 3 * t42;
t47 = t28 ^ 2;
t77 = 3 * t47;
t74 = pkin(2) + pkin(3);
t73 = pkin(3) - pkin(2);
t35 = pkin(2) ^ 2;
t12 = -t35 + t37;
t70 = t12 * t33;
t69 = t35 * t28;
t68 = t37 * t35;
t67 = t35 * qJ(1);
t66 = t41 * qJ(1);
t65 = pkin(2) - t20;
t64 = pkin(2) + t20;
t22 = -3 * t38;
t34 = t35 ^ 2;
t63 = -2 * t68 + 5 * t34 + t22;
t62 = t28 + t12;
t10 = pkin(1) - t65;
t11 = pkin(1) - t64;
t8 = pkin(1) + t64;
t9 = pkin(1) + t65;
t61 = t11 * t10 * t9 * t8;
t60 = (pkin(1) + t74) * (pkin(1) + t73) * (pkin(1) - t73) * (pkin(1) - t74);
t59 = -t28 - t35 - t82;
t58 = qJD(1) * (t37 + t59) / t8 ^ 2 / t9 ^ 2 / t10 ^ 2 / t11 ^ 2;
t57 = 1 / t20 * t80 * t58;
t30 = pkin(3) * t42;
t25 = qJ(1) * t47;
t23 = 3 * t37;
t1 = (((-15 * t42 + t63 + 18 * t70) * t28) + (t63 * t33) + ((-5 * t33 + t12) * t77) + (t12 * t78) + 0.6e1 * (-t25 + 0.2e1 * (-0.5e1 / 0.3e1 * t33 + t12) * t46 - (t42 - (2 * t70) + t38 + 0.2e1 / 0.3e1 * t68 - 0.5e1 / 0.3e1 * t34) * qJ(1)) * pkin(3) + t81) * t57;
t2 = [0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, -qJD(1) * (t59 * t79 - t64 * t65 * (-t42 - 4 * t66 + (-6 * t28 - 4 * t35) * t33 + (-4 * t46 - 8 * t67) * pkin(3) + t34 - 4 * t69 - t47) + (t38 + t78 + 12 * t66 + (18 * t28 - 2 * t35) * t33 + (12 * t46 - 4 * t67) * pkin(3) + 3 * t34 - 2 * t69 + t77) * t37) * (-t61) ^ (-0.1e1 / 0.2e1) * t80 / t61, 0, -((9 * t33 + 7 * t35 - 3 * t37) * t25 + (-5 * t42 + (46 * t35 - 6 * t37) * t33 + t79 + 6 * t68 - 9 * t34) * t46 - 9 * t30 * t28 + (5 * t47 + (34 * t35 + 6 * t37) * t28) * t41 + (5 * t49 + (29 * t35 - 9 * t37) * t47 + (-17 * t34 + t79 + 14 * t68) * t28 + (-t33 + t12) * t60) * pkin(3) + (t49 - (5 * t33 + t12) * t60) * qJ(1)) * t57, -(((t23 - 15 * t28 + t35) * t30) + ((-15 * t47 + (22 * t35 + 18 * t37) * t28 + t22 + 2 * t68 + t34) * t41) + ((4 * t35 * t62 - 20 * t42) * t46) + (-(6 * t44) + ((8 * t35 + 12 * t37) * t42) - 0.6e1 * (t47 + (-(2 * t37) - 0.14e2 / 0.3e1 * t35) * t28 + t38 - 0.4e1 / 0.3e1 * t68 + t34 / 0.3e1) * t33) * qJ(1) + (((t23 + 17 * t35) * t47 + (t22 - 7 * t34 + 10 * t68) * t28 + t81) * pkin(3))) * qJ(1) * t57, 0, 0, 0, 0, 0, -4 * t20 * (t62 + t82) * t58, 0, 0, 0, 0;];
cmat_reg = t2;
