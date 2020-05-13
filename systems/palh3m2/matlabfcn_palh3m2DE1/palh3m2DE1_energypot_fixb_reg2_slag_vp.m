% Calculate inertial parameters regressor of potential energy for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m2DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energypot_fixb_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:03:38
% EndTime: 2020-05-07 02:03:38
% DurationCPUTime: 0.27s
% Computational Cost: add. (215->64), mult. (408->99), div. (0->0), fcn. (427->22), ass. (0->54)
t229 = sin(qJ(1));
t235 = cos(qJ(1));
t206 = g(1) * t235 + g(2) * t229;
t234 = cos(qJ(2));
t211 = pkin(1) * t234 + pkin(12);
t228 = sin(qJ(2));
t243 = g(3) * (t228 * pkin(1) + pkin(11));
t249 = -t206 * t211 - t243;
t247 = g(3) * pkin(11);
t227 = sin(qJ(3));
t246 = pkin(4) * t227;
t242 = t227 * g(3);
t241 = t228 * t246 + pkin(12);
t233 = cos(qJ(3));
t212 = pkin(4) * t233 - pkin(1);
t240 = t212 * t228 - pkin(11);
t220 = sin(pkin(16));
t221 = cos(pkin(16));
t230 = sin(pkin(15));
t236 = cos(pkin(15));
t198 = t220 * t236 + t221 * t230;
t199 = -t220 * t230 + t221 * t236;
t188 = g(3) * t199 + t206 * t198;
t189 = -g(3) * t198 + t206 * t199;
t218 = pkin(17) + pkin(18);
t213 = sin(t218);
t214 = cos(t218);
t239 = t188 * t214 + t213 * t189;
t238 = t188 * t213 - t189 * t214;
t222 = sin(pkin(18));
t224 = cos(pkin(18));
t200 = t222 * t236 + t224 * t230;
t201 = -t222 * t230 + t224 * t236;
t237 = cos(pkin(14));
t232 = cos(qJ(4));
t231 = sin(pkin(14));
t226 = sin(qJ(4));
t225 = cos(pkin(17));
t223 = sin(pkin(17));
t219 = qJ(3) + qJ(2);
t216 = cos(t219);
t215 = sin(t219);
t207 = -g(1) * t229 + g(2) * t235;
t205 = t236 * pkin(8) - t230 * pkin(10);
t204 = t230 * pkin(8) + t236 * pkin(10);
t203 = t231 * t230 + t237 * t236;
t202 = -t237 * t230 + t231 * t236;
t197 = t206 * t233 + t242;
t196 = t233 * g(3) - t227 * t206;
t194 = -t220 * t204 + t205 * t221;
t193 = t204 * t221 + t220 * t205;
t191 = g(3) * t202 - t206 * t203;
t190 = g(3) * t203 + t206 * t202;
t1 = [0, 0, 0, 0, 0, 0, -t206, -t207, -g(3), -t247, 0, 0, 0, 0, 0, 0, -g(3) * t228 - t206 * t234, -g(3) * t234 + t206 * t228, t207, -t206 * pkin(12) - t247, 0, 0, 0, 0, 0, 0, t228 * t196 + t197 * t234, t196 * t234 - t228 * t197, t207, t249, 0, 0, 0, 0, 0, 0, -t238, t239, t207, (pkin(4) * t242 + t206 * t212) * t234 + g(3) * t240 - t206 * t241, 0, 0, 0, 0, 0, 0, t226 * t207 - t238 * t232, t232 * t207 + t238 * t226, -t239, (-t193 * g(3) + t206 * t194) * t214 + (-t194 * g(3) - t206 * t193) * t213 - g(3) * (-t234 * t246 - t240) - t206 * (-t212 * t234 + t241), 0, 0, 0, 0, 0, 0, -t190 * t228 + t191 * t234, -t190 * t234 - t228 * t191, t207, -g(3) * (pkin(11) + pkin(13)) + t206 * pkin(6), 0, 0, 0, 0, 0, 0, -g(3) * t200 + t206 * t201, g(3) * t201 + t206 * t200, t207, t249, 0, 0, 0, 0, 0, 0, g(3) * t215 + t206 * t216, g(3) * t216 - t206 * t215, t207, -t243 + t206 * ((-t200 * t223 + t201 * t225) * pkin(3) - t211) - g(3) * (t200 * t225 + t201 * t223) * pkin(3);];
U_reg = t1;
