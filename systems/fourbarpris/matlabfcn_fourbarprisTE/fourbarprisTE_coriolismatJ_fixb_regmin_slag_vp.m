% Calculate minimal parameter regressor of coriolis matrix for
% fourbarprisTE
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
% cmat_reg [(1*%NQJ)%x9]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbarprisTE_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:06
% EndTime: 2020-06-27 17:07:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (227->111), mult. (498->152), div. (24->7), fcn. (2->2), ass. (0->49)
t32 = (pkin(3) ^ 2);
t81 = 2 * qJ(1) * pkin(3) + t32;
t40 = pkin(3) * t32;
t43 = t40 ^ 2;
t27 = qJ(1) ^ 2;
t45 = qJ(1) * t27;
t48 = t45 ^ 2;
t75 = (pkin(1) - pkin(2));
t52 = (t75 ^ 2);
t74 = (pkin(1) + pkin(2));
t54 = (t74 ^ 2);
t80 = t75 * t52 * t74 * t54 - t43 - t48;
t19 = qJ(1) + pkin(3);
t79 = 1 / t19 ^ 2;
t36 = (pkin(1) ^ 2);
t37 = (t36 ^ 2);
t78 = 3 * t37;
t41 = t32 ^ 2;
t77 = 3 * t41;
t46 = t27 ^ 2;
t76 = 3 * t46;
t73 = pkin(2) + pkin(3);
t72 = pkin(3) - pkin(2);
t34 = pkin(2) ^ 2;
t11 = -t34 + t36;
t69 = t11 * t32;
t68 = t27 * t34;
t67 = t36 * t34;
t66 = t34 * qJ(1);
t65 = t40 * qJ(1);
t64 = pkin(2) - t19;
t63 = pkin(2) + t19;
t21 = -3 * t37;
t33 = t34 ^ 2;
t62 = -2 * t67 + 5 * t33 + t21;
t61 = t27 + t11;
t10 = pkin(1) - t63;
t7 = pkin(1) + t63;
t8 = pkin(1) + t64;
t9 = pkin(1) - t64;
t60 = t10 * t9 * t8 * t7;
t59 = (pkin(1) + t73) * (pkin(1) + t72) * (pkin(1) - t72) * (pkin(1) - t73);
t58 = -t27 - t34 - t81;
t57 = qJD(1) * (t36 + t58) / t7 ^ 2 / t8 ^ 2 / t9 ^ 2 / t10 ^ 2;
t56 = 1 / t19 * t79 * t57;
t29 = pkin(3) * t41;
t24 = qJ(1) * t46;
t22 = 3 * t36;
t1 = [(((-15 * t41 + t62 + 18 * t69) * t27) + (t62 * t32) + ((-5 * t32 + t11) * t76) + (t11 * t77) + 0.6e1 * (-t24 + 0.2e1 * (-0.5e1 / 0.3e1 * t32 + t11) * t45 - (t41 - (2 * t69) + t37 + 0.2e1 / 0.3e1 * t67 - 0.5e1 / 0.3e1 * t33) * qJ(1)) * pkin(3) + t80) * t56, 0, 0, -qJD(1) * (t58 * t78 - t63 * t64 * (-t41 - 4 * t65 + (-6 * t27 - 4 * t34) * t32 + (-4 * t45 - 8 * t66) * pkin(3) + t33 - 4 * t68 - t46) + (t37 + t77 + 12 * t65 + (18 * t27 - 2 * t34) * t32 + (12 * t45 - 4 * t66) * pkin(3) + 3 * t33 - 2 * t68 + t76) * t36) * (-t60) ^ (-0.1e1 / 0.2e1) * t79 / t60, -((9 * t32 + 7 * t34 - 3 * t36) * t24 + (-5 * t41 + (46 * t34 - 6 * t36) * t32 + t78 + 6 * t67 - 9 * t33) * t45 - 9 * t29 * t27 + (5 * t46 + (34 * t34 + 6 * t36) * t27) * t40 + (5 * t48 + (29 * t34 - 9 * t36) * t46 + (-17 * t33 + t78 + 14 * t67) * t27 + (-t32 + t11) * t59) * pkin(3) + (t48 - (5 * t32 + t11) * t59) * qJ(1)) * t56, -(((t22 - 15 * t27 + t34) * t29) + ((-15 * t46 + (22 * t34 + 18 * t36) * t27 + t21 + 2 * t67 + t33) * t40) + ((4 * t34 * t61 - 20 * t41) * t45) + (-(6 * t43) + ((8 * t34 + 12 * t36) * t41) - 0.6e1 * (t46 + (-(2 * t36) - 0.14e2 / 0.3e1 * t34) * t27 + t37 - 0.4e1 / 0.3e1 * t67 + t33 / 0.3e1) * t32) * qJ(1) + (((t22 + 17 * t34) * t46 + (t21 - 7 * t33 + 10 * t67) * t27 + t80) * pkin(3))) * qJ(1) * t56, -4 * t19 * (t61 + t81) * t57, 0, 0;];
cmat_reg = t1;
