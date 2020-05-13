% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = fourbar2OL_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:29
% EndTime: 2020-04-24 20:32:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (32->16), mult. (56->26), div. (0->0), fcn. (34->6), ass. (0->15)
t17 = -2 * pkin(2);
t11 = sin(qJ(1));
t14 = cos(qJ(1));
t16 = t14 * g(1) + t11 * g(2);
t15 = -t11 * g(1) + t14 * g(2);
t13 = cos(qJ(2));
t12 = cos(qJ(3));
t10 = sin(qJ(2));
t9 = sin(qJ(3));
t8 = qJDD(1) + qJDD(2);
t4 = -(qJD(1) ^ 2 * pkin(2)) - t16;
t3 = (qJDD(1) * pkin(2)) - t15;
t2 = (qJDD(1) + qJDD(2) / 0.2e1) * t17 + t15;
t1 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1) * t17 + t16;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t15, t16, 0, 0, 0, 0, 0, 0, 0, t8, -t1 * t10 + t2 * t13, -t1 * t13 - t10 * t2, 0, pkin(2) * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10 * t4 - t3 * t13, t10 * t3 + t13 * t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t9 * g(1) - t12 * g(2), t12 * g(1) + t9 * g(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauJ_reg = t5;
