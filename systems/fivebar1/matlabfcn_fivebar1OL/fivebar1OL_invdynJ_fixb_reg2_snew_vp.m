% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% fivebar1OL
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = fivebar1OL_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:17
% EndTime: 2020-04-27 06:13:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (60->34), mult. (112->56), div. (0->0), fcn. (60->8), ass. (0->33)
t17 = sin(qJ(4));
t36 = -0.2e1 * t17;
t19 = sin(qJ(2));
t35 = 0.2e1 * t19;
t34 = pkin(2) * (qJDD(1) + qJDD(2) / 0.2e1);
t33 = pkin(3) * (qJDD(3) + qJDD(4) / 0.2e1);
t20 = sin(qJ(1));
t32 = t20 * g(1);
t22 = cos(qJ(3));
t31 = t22 * g(2);
t12 = t22 * g(1);
t18 = sin(qJ(3));
t9 = t18 * g(2);
t30 = t12 + t9;
t11 = t20 * g(2);
t24 = cos(qJ(1));
t14 = t24 * g(1);
t29 = t14 + t11;
t28 = qJD(2) * pkin(2) * (qJD(1) + qJD(2) / 0.2e1);
t27 = qJD(4) * pkin(3) * (qJD(3) + qJD(4) / 0.2e1);
t13 = t24 * g(2);
t26 = -t13 + t32;
t10 = t18 * g(1);
t25 = t10 - t31;
t23 = cos(qJ(2));
t21 = cos(qJ(4));
t16 = qJDD(1) + qJDD(2);
t15 = qJDD(3) + qJDD(4);
t4 = qJD(3) ^ 2 * pkin(3) + t30;
t3 = -qJD(1) ^ 2 * pkin(2) - t29;
t2 = qJDD(1) * pkin(2) + t26;
t1 = qJDD(3) * pkin(3) + t25;
t5 = [0, 0, 0, 0, 0, qJDD(1), t26, t29, 0, 0, 0, 0, 0, 0, 0, t16, (-t26 - 0.2e1 * t34) * t23 + (-t14 / 0.2e1 - t11 / 0.2e1 + t28) * t35, (0.2e1 * t28 - t29) * t23 + (-t13 / 0.2e1 + t32 / 0.2e1 + t34) * t35, 0, pkin(2) * t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19 * t3 - t2 * t23, t19 * t2 + t23 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t25, t30, 0, 0, 0, 0, 0, 0, 0, t15, (t25 + 0.2e1 * t33) * t21 + (-t12 / 0.2e1 - t9 / 0.2e1 + t27) * t36, (-0.2e1 * t27 + t30) * t21 + (-t31 / 0.2e1 + t10 / 0.2e1 + t33) * t36, 0, pkin(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t21 * t1 + t17 * t4, -t17 * t1 + t4 * t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauJ_reg = t5;
