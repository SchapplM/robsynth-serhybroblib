% Jacobian of implicit kinematic constraints of
% fivebar1IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = fivebar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->20), mult. (38->29), div. (14->4), fcn. (22->11), ass. (0->15)
t42 = (qJ(3) + qJ(4));
t40 = qJ(1) - t42;
t31 = 0.1e1 / sin(qJ(2) + t40);
t39 = 0.1e1 / pkin(4);
t45 = t31 * t39;
t43 = qJ(1) + qJ(2);
t44 = cos((2 * t42)) - cos(0.2e1 * t43);
t38 = 0.1e1 / pkin(5);
t41 = t38 * t45;
t37 = 0.2e1 * qJ(1);
t36 = 2 * qJ(3);
t35 = sin(qJ(2));
t34 = sin(qJ(4));
t30 = 0.1e1 / t44;
t1 = [-(t44 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t36) + cos(qJ(2) + t37)) * pkin(2)) * t38 * t30, pkin(3) * t34 * t38 * t31; pkin(2) * t35 * t45, -(t44 * pkin(4) + (cos(qJ(4) + t36) - cos(qJ(4) + t37 + 0.2e1 * qJ(2))) * pkin(3)) * t39 * t30; pkin(2) * (-pkin(5) * t35 + pkin(4) * sin(t40)) * t41, pkin(3) * (pkin(4) * t34 + pkin(5) * sin(-qJ(3) + t43)) * t41;];
B21 = t1;
