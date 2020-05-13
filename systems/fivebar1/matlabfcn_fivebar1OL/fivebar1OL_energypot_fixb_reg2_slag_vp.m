% Calculate inertial parameters regressor of potential energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fivebar1OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:13
% EndTime: 2020-04-27 06:13:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->17), mult. (41->19), div. (0->0), fcn. (36->8), ass. (0->14)
t21 = pkin(1) * g(1);
t14 = sin(qJ(3));
t18 = cos(qJ(3));
t11 = g(1) * t18 + g(2) * t14;
t20 = cos(qJ(1));
t19 = cos(qJ(2));
t17 = cos(qJ(4));
t16 = sin(qJ(1));
t15 = sin(qJ(2));
t13 = sin(qJ(4));
t12 = g(1) * t20 + g(2) * t16;
t10 = g(1) * t16 - g(2) * t20;
t9 = g(1) * t14 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, -t12, t10, -g(3), 0, 0, 0, 0, 0, 0, 0, -t10 * t15 + t19 * t12, -t10 * t19 - t12 * t15, -g(3), -pkin(2) * t12, 0, 0, 0, 0, 0, 0, -t11, t9, -g(3), -t21, 0, 0, 0, 0, 0, 0, -t17 * t11 + t13 * t9, t13 * t11 + t9 * t17, -g(3), -t11 * pkin(3) - t21;];
U_reg = t1;
