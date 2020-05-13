% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = fivebar1OL_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:18
% EndTime: 2020-04-27 06:13:19
% DurationCPUTime: 0.42s
% Computational Cost: add. (163->69), mult. (198->52), div. (0->0), fcn. (128->8), ass. (0->37)
t223 = qJD(3) + qJD(4);
t219 = t223 ^ 2;
t221 = -qJDD(4) - qJDD(3);
t226 = sin(qJ(4));
t230 = cos(qJ(4));
t204 = t219 * t230 - t226 * t221;
t227 = sin(qJ(3));
t231 = cos(qJ(3));
t237 = t226 * t219 + t221 * t230;
t245 = t227 * t204 + t231 * t237;
t224 = qJD(1) + qJD(2);
t220 = t224 ^ 2;
t222 = qJDD(1) + qJDD(2);
t228 = sin(qJ(2));
t232 = cos(qJ(2));
t205 = t220 * t232 + t228 * t222;
t229 = sin(qJ(1));
t233 = cos(qJ(1));
t236 = t228 * t220 - t222 * t232;
t244 = t229 * t205 + t233 * t236;
t241 = -t233 * g(1) - t229 * g(2);
t240 = t227 * g(1) - t231 * g(2);
t239 = -t231 * g(1) - t227 * g(2);
t238 = t229 * g(1) - t233 * g(2);
t235 = qJD(1) ^ 2;
t213 = t233 * qJDD(1) - t229 * t235;
t211 = -t229 * qJDD(1) - t233 * t235;
t234 = qJD(3) ^ 2;
t212 = t231 * qJDD(3) - t227 * t234;
t210 = -t227 * qJDD(3) - t231 * t234;
t209 = -t235 * pkin(2) + t241;
t208 = -t234 * pkin(3) + t239;
t207 = qJDD(1) * pkin(2) + t238;
t206 = qJDD(3) * pkin(3) + t240;
t199 = t205 * t233 - t229 * t236;
t198 = -t204 * t231 + t227 * t237;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t211, -t213, 0, -g(1), 0, 0, 0, 0, 0, 0, t199, -t244, 0, t211 * pkin(2) - g(1), 0, 0, 0, 0, 0, 0, t210, -t212, 0, -g(1), 0, 0, 0, 0, 0, 0, t198, t245, 0, t210 * pkin(3) - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t213, t211, 0, -g(2), 0, 0, 0, 0, 0, 0, t244, t199, 0, t213 * pkin(2) - g(2), 0, 0, 0, 0, 0, 0, t212, t210, 0, -g(2), 0, 0, 0, 0, 0, 0, -t245, t198, 0, t212 * pkin(3) - g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, -qJDD(1), 0, t241, 0, 0, 0, 0, 0, 0, t205, -t236, 0, t209, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t235, 0, t238, 0, 0, 0, 0, 0, 0, t236, t205, 0, t207, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, -t222, 0, -t228 * t207 - t209 * t232, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t220, 0, -t207 * t232 + t228 * t209, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -qJDD(3), 0, t239, 0, 0, 0, 0, 0, 0, -t204, t237, 0, t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t234, 0, t240, 0, 0, 0, 0, 0, 0, -t237, -t204, 0, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t221, 0, t226 * t206 + t230 * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, -t219, 0, t230 * t206 - t226 * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
f_new_reg = t1;
