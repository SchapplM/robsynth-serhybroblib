% Calculate minimal parameter regressor of potential energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh2m2OL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:35
% EndTime: 2020-06-30 18:06:36
% DurationCPUTime: 0.05s
% Computational Cost: add. (80->22), mult. (64->34), div. (0->0), fcn. (68->12), ass. (0->22)
t165 = qJ(2) + qJ(3);
t164 = qJ(4) + t165;
t161 = qJ(5) + t164;
t157 = sin(t161);
t177 = g(3) * t157;
t166 = sin(qJ(6));
t168 = sin(qJ(1));
t176 = t168 * t166;
t169 = cos(qJ(6));
t175 = t168 * t169;
t171 = cos(qJ(1));
t174 = t171 * t166;
t173 = t171 * t169;
t172 = g(1) * t171 + g(2) * t168;
t170 = cos(qJ(2));
t167 = sin(qJ(2));
t163 = cos(t165);
t162 = sin(t165);
t160 = cos(t164);
t159 = sin(t164);
t158 = cos(t161);
t1 = [0, -t172, g(1) * t168 - g(2) * t171, 0, 0, 0, 0, 0, -g(3) * t167 - t172 * t170, -g(3) * t170 + t172 * t167, 0, 0, 0, 0, 0, -g(3) * t162 - t172 * t163, -g(3) * t163 + t172 * t162, 0, 0, 0, 0, 0, -g(3) * t159 - t172 * t160, -g(3) * t160 + t172 * t159, 0, 0, 0, 0, 0, -t172 * t158 - t177, -g(3) * t158 + t172 * t157, 0, 0, 0, 0, 0, -g(1) * (t158 * t173 - t176) - g(2) * (t158 * t175 + t174) - t169 * t177, -g(1) * (-t158 * t174 - t175) - g(2) * (-t158 * t176 + t173) + t166 * t177;];
U_reg = t1;
