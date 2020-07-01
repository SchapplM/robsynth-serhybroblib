% Calculate minimal parameter regressor of potential energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m1OL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:10
% EndTime: 2020-06-30 17:45:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (52->19), mult. (54->30), div. (0->0), fcn. (58->10), ass. (0->19)
t135 = qJ(2) + qJ(3);
t134 = qJ(4) + t135;
t130 = sin(t134);
t147 = g(3) * t130;
t136 = sin(qJ(5));
t138 = sin(qJ(1));
t146 = t138 * t136;
t139 = cos(qJ(5));
t145 = t138 * t139;
t141 = cos(qJ(1));
t144 = t141 * t136;
t143 = t141 * t139;
t142 = g(1) * t141 + g(2) * t138;
t140 = cos(qJ(2));
t137 = sin(qJ(2));
t133 = cos(t135);
t132 = sin(t135);
t131 = cos(t134);
t1 = [0, -t142, g(1) * t138 - g(2) * t141, 0, 0, 0, 0, 0, g(3) * t137 - t142 * t140, g(3) * t140 + t142 * t137, 0, 0, 0, 0, 0, g(3) * t132 - t142 * t133, g(3) * t133 + t142 * t132, 0, 0, 0, 0, 0, -t142 * t131 + t147, g(3) * t131 + t142 * t130, 0, 0, 0, 0, 0, -g(1) * (t131 * t143 - t146) - g(2) * (t131 * t145 + t144) + t139 * t147, -g(1) * (-t131 * t144 - t145) - g(2) * (-t131 * t146 + t143) - t136 * t147;];
U_reg = t1;
