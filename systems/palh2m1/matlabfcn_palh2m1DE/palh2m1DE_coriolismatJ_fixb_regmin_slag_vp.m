% Calculate minimal parameter regressor of coriolis matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = palh2m1DE_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:38
% EndTime: 2020-06-30 17:39:38
% DurationCPUTime: 0.30s
% Computational Cost: add. (236->74), mult. (495->116), div. (0->0), fcn. (371->6), ass. (0->79)
t44 = cos(qJ(2));
t38 = t44 ^ 2;
t87 = 0.2e1 * t38;
t43 = cos(qJ(3));
t73 = t44 * t43;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t76 = t41 * t40;
t24 = t73 - t76;
t35 = qJD(2) + qJD(3);
t21 = t35 * t24;
t74 = t44 * t40;
t75 = t43 * t41;
t23 = t74 + t75;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t65 = qJD(1) * t42;
t86 = (-t39 * t21 + t23 * t65) * pkin(3);
t66 = qJD(1) * t39;
t85 = (-t42 * t21 - t23 * t66) * pkin(3);
t84 = t87 - 0.1e1;
t83 = pkin(2) * (qJD(2) + qJD(3) / 0.2e1);
t82 = pkin(3) * t23;
t81 = t44 * pkin(2);
t80 = t35 * t40;
t79 = t35 * t43;
t36 = qJD(1) + qJD(4);
t78 = t36 * t39;
t77 = t36 * t42;
t72 = pkin(2) * qJD(2);
t71 = pkin(2) * qJD(3);
t37 = t43 ^ 2;
t9 = t41 * (t37 - 0.1e1 / 0.2e1) * t44 + (t38 - 0.1e1 / 0.2e1) * t43 * t40;
t70 = qJD(1) * t9;
t26 = pkin(1) * t44 + pkin(2) * t87 - pkin(2);
t30 = pkin(1) + 0.2e1 * t81;
t69 = qJD(1) * (t26 * t43 - t30 * t76);
t16 = -0.4e1 * t73 * t76 + (0.4e1 * t38 - 0.2e1) * t37 - t84;
t68 = qJD(1) * t16;
t67 = qJD(1) * (pkin(1) + t81);
t64 = qJD(1) * t44;
t63 = qJD(2) * t41;
t62 = qJD(2) * t44;
t29 = t43 * pkin(3) + pkin(2);
t19 = -pkin(3) * t76 + t29 * t44 + pkin(1) + pkin(4);
t61 = qJD(4) * t19;
t60 = (t26 * t40 + t30 * t75) * qJD(1);
t55 = pkin(3) * t74;
t22 = t41 * t29 + t55;
t59 = t22 * qJD(1);
t58 = pkin(1) * t79;
t57 = t40 * t83;
t56 = t36 * t82;
t54 = t22 * t78;
t53 = t35 * t76;
t52 = pkin(1) * t41 * qJD(1);
t51 = pkin(1) * t64;
t50 = qJD(1) * t82;
t49 = t23 * t67;
t48 = t24 * t67;
t47 = t41 * t64;
t46 = -0.2e1 * t43 * t83;
t45 = t39 * t56;
t25 = pkin(3) * t79 + t72;
t12 = t41 * t25 + t35 * t55;
t33 = t43 * t72;
t32 = t40 * t72;
t27 = t84 * qJD(1);
t20 = t23 * t35;
t18 = t42 * t56;
t17 = t22 * t77;
t15 = -pkin(3) * t53 + t25 * t44;
t8 = t12 * t42;
t7 = t39 * t12;
t6 = -0.2e1 * t70;
t5 = 0.2e1 * t70;
t2 = t15 * t42 + t39 * t59;
t1 = t15 * t39 - t42 * t59;
t3 = [0, 0, 0, t41 * t62, t84 * qJD(2), 0, 0, 0, -pkin(1) * t63, -pkin(1) * t62, 0.2e1 * t35 * t9, t35 * t16, 0, 0, 0, -0.2e1 * t38 * t57 + (-pkin(1) * t80 + t41 * t46) * t44 + t32 - t41 * t58, t38 * t46 + (0.2e1 * t41 * t57 - t58) * t44 + t33 + pkin(1) * t53, -t12, 0, -t39 * t61 - t8, -t42 * t61 + t7; 0, 0, 0, t47, t27, -t62, t63, 0, -t52, -t51, t5, t68, -t21, t20, 0, -t60, -t69, -t59, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t68, -t21, t20, 0, -t49, -t48, -t50, 0, -t86, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t78, -t19 * t77; 0, 0, 0, -t47, -t27, 0, 0, 0, t52, t51, t6, -t68, 0, 0, 0, t60, t69, t59, 0, t17, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40 * t71, -t43 * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(2) * t80, -pkin(2) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t68, 0, 0, 0, t49, t48, t50, 0, t18, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t66 - t8, t19 * t65 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
