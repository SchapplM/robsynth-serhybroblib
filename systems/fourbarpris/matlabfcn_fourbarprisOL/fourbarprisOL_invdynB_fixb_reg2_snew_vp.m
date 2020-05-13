% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% fourbarprisOL
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% tauB_reg [6x(4*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = fourbarprisOL_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynB_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:09
% EndTime: 2020-05-07 09:52:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (122->60), mult. (239->48), div. (0->0), fcn. (168->4), ass. (0->33)
t55 = qJD(1) ^ 2;
t54 = qJD(3) ^ 2;
t41 = cos(qJ(3));
t53 = t41 * g(3);
t42 = cos(qJ(1));
t52 = t42 * g(3);
t51 = qJ(2) * g(3);
t50 = qJDD(1) * qJ(2);
t40 = sin(qJ(1));
t25 = t40 * g(1) - t42 * g(2);
t49 = (2 * qJD(2) * qJD(1)) - t25;
t27 = t42 * g(1) + t40 * g(2);
t17 = t49 + t50;
t18 = -t55 * qJ(2) + qJDD(2) + t27;
t15 = -t42 * t17 - t40 * t18;
t48 = t40 * t17 - t42 * t18;
t39 = sin(qJ(3));
t24 = t39 * g(1) - t41 * g(2);
t26 = t41 * g(1) + t39 * g(2);
t47 = t41 * t24 - t39 * t26;
t46 = -t39 * t24 - t41 * t26;
t16 = t42 * t25 - t40 * t27;
t45 = -t40 * t25 - t42 * t27;
t22 = -t40 * qJDD(1) - t42 * t55;
t44 = -pkin(1) * t22 - t27;
t43 = pkin(1) * g(3);
t35 = t40 * g(3);
t34 = t39 * g(3);
t23 = t42 * qJDD(1) - t40 * t55;
t21 = -t41 * qJDD(3) + t39 * t54;
t20 = -t39 * qJDD(3) - t41 * t54;
t19 = pkin(1) * t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t22, t23, 0, t45, 0, 0, 0, 0, 0, 0, -t23, 0, -t22, t48, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, -t23, -t22, 0, t16, 0, 0, 0, 0, 0, 0, t22, 0, -t23, t15, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, -t23, 0, -t22, 0, t35, t52, -t16, 0, 0, -t22, 0, 0, t23, 0, -t52, t15, t35, t40 * t51, 0, 0, t21, 0, -t20, 0, t34, t53, -t47, 0; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t22, 0, -t23, 0, -t52, t35, t45, t43, 0, -t23, 0, 0, -t22, 0, -t35, -t48, -t52, -t42 * t51 + t43, 0, 0, t20, 0, t21, 0, -t53, t34, t46, 0; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t19 - t25, t44, 0, pkin(1) * t16, 0, 0, 0, qJDD(1), 0, 0, qJDD(2) - t44, 0, -t19 + t49 + 0.2e1 * t50, pkin(1) * t15 + qJ(2) * t17, 0, 0, 0, 0, 0, qJDD(3), -t24, -t26, 0, 0;];
tauB_reg = t1;
