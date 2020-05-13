% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
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
% tauB_reg [6x(4*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = fourbar2OL_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynB_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:29
% EndTime: 2020-04-24 20:32:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (133->65), mult. (202->57), div. (0->0), fcn. (136->6), ass. (0->35)
t54 = -2 * pkin(2);
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t42 = cos(qJ(2));
t43 = cos(qJ(1));
t53 = (-t40 * t39 + t42 * t43) * g(3);
t38 = sin(qJ(3));
t52 = t38 * g(3);
t51 = t40 * g(3);
t41 = cos(qJ(3));
t50 = t41 * g(3);
t49 = t43 * g(3);
t48 = -t40 * g(1) + t43 * g(2);
t47 = t43 * g(1) + t40 * g(2);
t37 = qJD(1) + qJD(2);
t35 = t37 ^ 2;
t36 = qJDD(1) + qJDD(2);
t22 = t39 * t35 - t36 * t42;
t23 = t35 * t42 + t39 * t36;
t18 = t22 * t43 + t23 * t40;
t46 = t22 * t40 - t23 * t43;
t45 = qJD(1) ^ 2;
t31 = t43 * qJDD(1) - t40 * t45;
t29 = t40 * qJDD(1) + t43 * t45;
t44 = qJD(3) ^ 2;
t33 = -pkin(1) * t44 + g(1);
t32 = pkin(1) * qJDD(3) - g(2);
t30 = t41 * qJDD(3) - t38 * t44;
t28 = t38 * qJDD(3) + t41 * t44;
t26 = t31 * pkin(2) - g(2);
t25 = -t29 * pkin(2) - g(1);
t24 = (t39 * t43 + t40 * t42) * g(3);
t21 = (qJDD(1) + qJDD(2) / 0.2e1) * t54 + t48;
t20 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1) * t54 + t47;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t29, -t31, 0, -g(1), 0, 0, 0, 0, 0, 0, -t46, -t18, 0, t25, 0, 0, 0, 0, 0, 0, -t28, -t30, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t31, -t29, 0, -g(2), 0, 0, 0, 0, 0, 0, t18, -t46, 0, t26, 0, 0, 0, 0, 0, 0, t30, -t28, 0, -g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t31, 0, -t29, 0, -t51, -t49, g(2), 0, 0, 0, t18, 0, -t46, 0, t24, t53, -t26, -pkin(2) * t51, 0, 0, t30, 0, -t28, 0, -t52, -t50, g(2), 0; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t29, 0, t31, 0, t49, -t51, -g(1), 0, 0, 0, t46, 0, t18, 0, -t53, t24, t25, pkin(2) * t49, 0, 0, t28, 0, t30, 0, t50, -t52, -g(1), pkin(1) * g(3); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t48, t47, 0, 0, 0, 0, 0, 0, 0, t36, -t20 * t39 + t21 * t42, -t20 * t42 - t39 * t21, 0, pkin(2) * (qJDD(1) * pkin(2) - t48), 0, 0, 0, 0, 0, qJDD(3), t32 * t41 + t33 * t38, -t32 * t38 + t33 * t41, 0, -pkin(1) * g(2);];
tauB_reg = t1;
