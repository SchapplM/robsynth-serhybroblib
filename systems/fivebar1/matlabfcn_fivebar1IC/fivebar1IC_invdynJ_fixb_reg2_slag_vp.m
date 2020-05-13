% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% tau_reg [2x(5*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fivebar1IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:10
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (196->56), mult. (303->76), div. (24->4), fcn. (172->15), ass. (0->46)
t30 = sin(qJ(3));
t34 = cos(qJ(3));
t12 = g(1) * t30 - g(2) * t34;
t61 = qJDD(3) * pkin(3) + t12;
t32 = sin(qJ(1));
t58 = cos(qJ(1));
t13 = g(1) * t32 - g(2) * t58;
t60 = qJDD(1) * pkin(2) + t13;
t29 = sin(qJ(4));
t59 = pkin(3) * t29;
t51 = qJ(4) + qJ(3);
t52 = qJ(1) + qJ(2);
t57 = cos(0.2e1 * t51) - cos(0.2e1 * t52);
t28 = qJD(1) + qJD(2);
t56 = qJD(1) * t28;
t55 = qJD(2) * t28;
t36 = 0.2e1 * qJ(3);
t37 = 0.2e1 * qJ(1);
t38 = 0.1e1 / pkin(5);
t7 = 0.1e1 / t57;
t50 = t38 * (t57 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t36) + cos(qJ(2) + t37)) * pkin(2)) * t7;
t39 = 0.1e1 / pkin(4);
t49 = t39 * (t57 * pkin(4) + (cos(qJ(4) + t36) - cos(qJ(4) + t37 + 0.2e1 * qJ(2))) * pkin(3)) * t7;
t26 = qJDD(1) + qJDD(2);
t25 = qJDD(3) + qJDD(4);
t31 = sin(qJ(2));
t35 = cos(qJ(2));
t47 = pkin(2) * qJD(1) * qJD(2);
t48 = t60 * t31 + t35 * t47;
t14 = g(1) * t34 + g(2) * t30;
t33 = cos(qJ(4));
t46 = -t29 * t12 + t33 * t14;
t16 = -0.1e1 / sin(t51 - t52);
t45 = pkin(2) * t16 * t31 * t39;
t44 = t16 * t38 * t59;
t27 = qJD(3) + qJD(4);
t43 = qJD(3) * (-qJD(4) + t27);
t42 = qJD(4) * (-qJD(3) - t27);
t41 = t29 * t14 + t61 * t33;
t15 = g(1) * t58 + g(2) * t32;
t40 = -t35 * t13 + (-t15 + t47) * t31;
t4 = (-pkin(2) * t56 - t15) * t35 + t48;
t3 = (-qJDD(1) * t35 - t31 * t56) * pkin(2) + t40;
t2 = (-qJDD(3) * t29 + t33 * t43) * pkin(3) + t46;
t1 = t43 * t59 + t41;
t5 = [0, 0, 0, 0, 0, qJDD(1), t13, t15, 0, 0, 0, 0, 0, 0, 0, (-t50 + 0.1e1) * t26, -t3 * t50 + (t31 * t55 + (-qJDD(1) - t26) * t35) * pkin(2) + t40, -t4 * t50 - t35 * t15 + (t26 * t31 + t35 * t55) * pkin(2) + t48, 0, t60 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t45, t1 * t45, t2 * t45, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t44, t3 * t44, t4 * t44, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t12, t14, 0, 0, 0, 0, 0, 0, 0, (-t49 + 0.1e1) * t25, -t1 * t49 + (t25 * t33 + t29 * t42) * pkin(3) + t41, -t2 * t49 + ((-qJDD(3) - t25) * t29 + t33 * t42) * pkin(3) + t46, 0, t61 * pkin(3);];
tau_reg = t5;
